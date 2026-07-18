import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import os
import io
import math

st.set_page_config(page_title="AS-MS Pooling Engine", page_icon="🧪", layout="wide")

st.title("NCATS ASMS Advanced Pooling Engine")
st.markdown("Ingests structural `.sdf` library exports to predict ionization modes, optimize mass-diversity using a modular stride algorithm, calculate volume-normalizing DMSO back-flushes, and determine exact Echo nanoliter dispense steps.")

# ==========================================
# 1. Sidebar Control Panel
# ==========================================
st.sidebar.header("Configuration Panel")

st.sidebar.subheader("1. Library Pooling Options")
pool_size = st.sidebar.number_input("Target Compounds per Well", min_value=2, max_value=50, value=10, step=1)
min_mz_threshold = st.sidebar.number_input("Minimum Allowed m/z Delta (Da)", min_value=0.5, max_value=10.0, value=2.0, step=0.5)

st.sidebar.subheader("2. Volumetric Normalization")
vol_per_comp = st.sidebar.number_input("Source Plate: Vol per Compound (nL)", min_value=10, max_value=2000, value=100, step=50, help="The nanoliter volume of each individual compound added to create a pooled source well.")

st.sidebar.subheader("3. Echo Assay Calculator")
dest_well_vol = st.sidebar.number_input("Assay Plate: Total Well Volume (µL)", min_value=1.0, max_value=500.0, value=50.0, step=5.0)
desired_conc = st.sidebar.number_input("Assay Plate: Target Concentration (µM)", min_value=0.1, max_value=100.0, value=10.0, step=1.0)

# File Uploader
st.markdown("### Upload Library File")
uploaded_file = st.file_uploader("Choose an SDF file", type=["sdf"])
plate_prefix = st.text_input("Plate Name Prefix", value="ASMS")

# ==========================================
# 2. Computational Core Functions
# ==========================================
def process_sdf(file_path):
    supplier = Chem.SDMolSupplier(file_path)
    data = []
    
    basic_nitrogen = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(N-[#6a])]")
    acidic_group = Chem.MolFromSmarts("[C,S](=[O,S])[O;H1,-1]")
    
    for idx, mol in enumerate(supplier):
        if mol is None: continue
        
        sample_id = None
        for prop_name in ['SAMPLE_ID', 'Name', 'ID', 'sample_id', 'id']:
            if mol.HasProp(prop_name):
                sample_id = mol.GetProp(prop_name)
                break
        
        if not sample_id: 
            sample_id = f"UNKNOWN_{idx}"
        else:
            if '-' in str(sample_id):
                sample_id = str(sample_id).split('-')[0]
            
        try:
            exact_mass = Descriptors.ExactMolWt(mol)
            smiles = Chem.MolToSmiles(mol)
            
            if exact_mass == 0 or not smiles or smiles.strip() == "":
                continue
                
            has_base = mol.HasSubstructMatch(basic_nitrogen)
            has_acid = mol.HasSubstructMatch(acidic_group)
            
            if has_base and not has_acid:
                predicted_mode = "positive"
                target_mz = exact_mass + 1.0073
            elif has_acid and not has_base:
                predicted_mode = "negative"
                target_mz = exact_mass - 1.0073
            else:
                predicted_mode = "positive"  # Amphoteric/Neutral default
                target_mz = exact_mass + 1.0073
                
            data.append({
                'SAMPLE_ID': sample_id,
                'Exact_Mass': exact_mass,
                'Target_m_z': target_mz,
                'Predicted_Mode': predicted_mode,
                'SMILES': smiles
            })
        except:
            continue
    return pd.DataFrame(data)

def assign_wells_advanced(df, target_size, prefix, vol_comp, assay_vol, assay_conc):
    pooled_records = []
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    columns = range(1, 25)
    well_coordinates = [f"{r}{c:02d}" for r in rows for c in columns]
    
    current_plate = 1
    well_pointer = 0
    clean_prefix = prefix.strip().rstrip('_')
    
    # Calculate Echo Transfer Step: Vol = (Conc * Dest_Vol) / 10
    echo_transfer_nl = (assay_conc * assay_vol) / 10.0
    
    for mode, group in df.groupby('Predicted_Mode'):
        sorted_group = group.sort_values(by='Target_m_z').reset_index(drop=True)
        num_compounds = len(sorted_group)
        
        if num_compounds == 0: continue
        
        # Calculate optimal number of pools using ceiling division to never drop remainders
        num_pools = math.ceil(num_compounds / target_size)
        stride_pools = [[] for _ in range(num_pools)]
        
        for idx, row in sorted_group.iterrows():
            pool_idx = idx % num_pools
            stride_pools[pool_idx].append(row)
            
        for pool in stride_pools:
            if not pool: continue
            if well_pointer >= len(well_coordinates):
                well_pointer = 0
                current_plate += 1
                
            assigned_well = well_coordinates[well_pointer]
            well_pointer += 1
            
            actual_pool_count = len(pool)
            # Volume calculations for fluidic normalization
            dmso_backflush_nl = (target_size - actual_pool_count) * vol_comp
            total_well_fluid_volume_nl = target_size * vol_comp
            
            # Calculate actual minimum m/z delta within this specific well pool
            pool_mzs = sorted([comp['Target_m_z'] for comp in pool])
            if len(pool_mzs) > 1:
                min_delta_observed = min(pool_mzs[i+1] - pool_mzs[i] for i in range(len(pool_mzs)-1))
            else:
                min_delta_observed = float('inf')
                
            for sub_idx, comp in enumerate(pool):
                pooled_records.append({
                    'Destination_Plate': f"{clean_prefix}_PLT_{current_plate}",
                    'Destination_Well': assigned_well,
                    'Well_Sub_Index': sub_idx + 1,
                    'Compounds_In_Pool': actual_pool_count,
                    'Backflush_Required': "YES" if dmso_backflush_nl > 0 else "NO",
                    'NCGC_ID': comp['SAMPLE_ID'],
                    'Exact_Mass': round(comp['Exact_Mass'], 4),
                    'Target_m_z': round(comp['Target_m_z'], 4),
                    'Min_m_z_Delta_In_Well': round(min_delta_observed, 4) if min_delta_observed != float('inf') else 0.0,
                    'DMSO_Backflush_Volume_nL': dmso_backflush_nl,
                    'Total_Well_Fluid_Vol_nL': total_well_fluid_volume_nl,
                    'Echo_Transfer_To_Destination_nL': echo_transfer_nl,
                    'Ionization_Mode': comp['Predicted_Mode'],
                    'SMILES': comp['SMILES']
                })
                
    return pd.DataFrame(pooled_records)

def generate_interactive_html(df):
    import json
    plate_data_dict = {}
    for _, row in df.iterrows():
        plt = row['Destination_Plate']
        well = row['Destination_Well']
        
        if plt not in plate_data_dict: plate_data_dict[plt] = {}
        if well not in plate_data_dict[plt]: plate_data_dict[plt][well] = []
            
        svg_text = ""
        try:
            mol = Chem.MolFromSmiles(row['SMILES'])
            if mol:
                drawer = rdMolDraw2D.MolDraw2DSVG(160, 160)
                clean_mol = rdMolDraw2D.PrepareMolForDrawing(mol)
                drawer.DrawMolecule(clean_mol)
                drawer.FinishDrawing()
                svg_text = drawer.GetDrawingText()
        except:
            svg_text = ""
            
        plate_data_dict[plt][well].append({
            'id': row['NCGC_ID'],
            'mass': row['Exact_Mass'],
            'mz': row['Target_m_z'],
            'smiles': row['SMILES'],
            'img': svg_text,
            'mode': row['Ionization_Mode'],
            'backflush': int(row['DMSO_Backflush_Volume_nL'])
        })
        
    js_data_payload = json.dumps(plate_data_dict)
    
    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>NCATS ASMS Interactive Plate Mapper</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; background-color: #f8f9fa; margin: 0; padding: 20px; color: #333; }}
        .header {{ display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid #e9ecef; padding-bottom: 15px; margin-bottom: 20px; }}
        h1 {{ margin: 0; font-size: 24px; color: #1e293b; }}
        select {{ padding: 8px 16px; font-size: 16px; border-radius: 6px; border: 1px solid #cbd5e1; background: white; cursor: pointer; }}
        .main-container {{ display: flex; gap: 24px; align-items: flex-start; }}
        .plate-box {{ background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); border: 1px solid #e2e8f0; }}
        .map-legend {{ display: flex; gap: 24px; margin-bottom: 20px; font-size: 13px; font-weight: 600; justify-content: center; background: #f8fafc; padding: 10px; border-radius: 8px; border: 1px solid #e2e8f0; flex-wrap: wrap; }}
        .legend-item {{ display: flex; align-items: center; gap: 8px; }}
        .grid-384 {{ display: grid; grid-template-columns: 30px repeat(24, 26px); gap: 4px; align-items: center; justify-items: center; }}
        .col-header {{ font-size: 11px; font-weight: bold; color: #64748b; text-align: center; width: 26px; }}
        .row-header {{ font-size: 11px; font-weight: bold; color: #64748b; text-align: center; height: 26px; display: flex; align-items: center; justify-content: center; }}
        .well {{ width: 22px; height: 22px; border-radius: 50%; background-color: #f1f5f9; border: 1px solid #cbd5e1; cursor: pointer; transition: all 0.15s ease; box-sizing: border-box; }}
        .well.populated.positive {{ background-color: #bfdbfe; border-color: #3b82f6; }}
        .well.populated.negative {{ background-color: #fecdd3; border-color: #f43f5e; }}
        .well.backflush-needed {{ border-style: dashed !important; border-width: 2px !important; }}
        .well:hover {{ transform: scale(1.2); border-color: #475569 !important; box-shadow: 0 0 4px rgba(0,0,0,0.2); z-index: 10; }}
        .well.active-well {{ border-color: #1e3a8a !important; background-color: #eff6ff !important; box-shadow: 0 0 0 3px #3b82f6; }}
        .display-panel {{ flex: 1; min-width: 400px; background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); border: 1px solid #e2e8f0; max-height: 85vh; overflow-y: auto; }}
        .panel-title {{ font-size: 18px; font-weight: bold; margin-bottom: 15px; color: #1e293b; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 10px; }}
        .backflush-tag {{ font-size: 12px; font-weight: bold; color: #b45309; background-color: #fef3c7; padding: 4px 10px; border-radius: 6px; border: 1px solid #fde68a; }}
        .compound-card {{ display: flex; align-items: center; gap: 15px; padding: 12px; margin-bottom: 12px; border: 1px solid #e2e8f0; border-radius: 8px; background: #f8fafc; }}
        .compound-info {{ flex: 1; font-size: 13px; }}
        .compound-id {{ font-size: 15px; font-weight: bold; color: #2563eb; margin-bottom: 4px; }}
        .struct-img {{ width: 130px; height: 130px; background: white; border: 1px solid #e2e8f0; border-radius: 6px; display: flex; align-items: center; justify-content: center; overflow: hidden; }}
        .struct-img svg {{ max-width: 100%; max-height: 100%; }}
        .placeholder-text {{ color: #94a3b8; font-style: italic; text-align: center; margin-top: 50px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>NCATS ASMS Interactive Pool Navigator</h1>
        <select id="plateSelect" onchange="renderPlate()"></select>
    </div>
    <div class="main-container">
        <div class="plate-box">
            <div class="map-legend">
                <div class="legend-item">
                    <div style="width: 14px; height: 14px; border-radius: 50%; background-color: #fecdd3; border: 1px solid #f43f5e;"></div>
                    <span>Negative Pools</span>
                </div>
                <div class="legend-item">
                    <div style="width: 14px; height: 14px; border-radius: 50%; background-color: #bfdbfe; border: 1px solid #3b82f6;"></div>
                    <span>Positive Pools</span>
                </div>
                <div class="legend-item">
                    <div style="width: 14px; height: 14px; border-radius: 50%; background-color: #f1f5f9; border: 2px dashed #475569;"></div>
                    <span>Dashed Border = Requires DMSO Back-flush</span>
                </div>
            </div>
            <div id="gridContainer" class="grid-384"></div>
        </div>
        <div class="display-panel">
            <div id="panelTitle" class="panel-title">Well Pool Inspector</div>
            <div id="compoundsContainer">
                <div class="placeholder-text">Click any populated well on the left layout grid to inspect its contents...</div>
            </div>
        </div>
    </div>
    <script>
        const db = {js_data_payload};
        const rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'];
        const select = document.getElementById('plateSelect');
        Object.keys(db).sort().forEach(plt => {{
            let opt = document.createElement('option');
            opt.value = plt; opt.innerHTML = plt; select.appendChild(opt);
        }});
        function renderPlate() {{
            const currentPlate = select.value;
            const container = document.getElementById('gridContainer');
            container.innerHTML = '';
            let spacer = document.createElement('div'); container.appendChild(spacer);
            for(let c=1; c<=24; c++) {{
                let header = document.createElement('div'); header.className = 'col-header';
                header.innerHTML = c < 10 ? '0'+c : c; container.appendChild(header);
            }}
            rows.forEach(r => {{
                let rHeader = document.createElement('div'); rHeader.className = 'row-header';
                rHeader.innerHTML = r; container.appendChild(rHeader);
                for(let c=1; c<=24; c++) {{
                    let wellName = r + (c < 10 ? '0'+c : c);
                    let wellDiv = document.createElement('div'); wellDiv.className = 'well'; wellDiv.id = wellName;
                    const dynamicData = db[currentPlate] && db[currentPlate][wellName];
                    if (dynamicData) {{
                        wellDiv.classList.add('populated');
                        const wellMode = dynamicData[0].mode;
                        wellDiv.classList.add(wellMode);
                        if(dynamicData[0].backflush > 0) {{
                            wellDiv.classList.add('backflush-needed');
                        }}
                        wellDiv.onclick = () => selectWell(wellName, dynamicData, wellDiv);
                    }}
                    container.appendChild(wellDiv);
                }}
            }});
        }}
        function selectWell(wellName, compounds, element) {{
            document.querySelectorAll('.well').forEach(w => w.classList.remove('active-well'));
            element.classList.add('active-well');
            
            const wellModeLabel = compounds[0].mode.toUpperCase();
            let headerHTML = `<span>Contents of Well: ${{wellName}} (${{wellModeLabel}} Mode)</span>`;
            if(compounds[0].backflush > 0) {{
                headerHTML += `<span class="backflush-tag">⚠️ DMSO Back-flush: +${{compounds[0].backflush}} nL</span>`;
            }}
            document.getElementById('panelTitle').innerHTML = headerHTML;
            
            const listContainer = document.getElementById('compoundsContainer');
            listContainer.innerHTML = '';
            compounds.forEach(c => {{
                let card = document.createElement('div'); card.className = 'compound-card';
                let info = document.createElement('div'); info.className = 'compound-info';
                info.innerHTML = `
                    <div class="compound-id">${{c.id}}</div>
                    <div><strong>Exact Mass:</strong> ${{c.mass.toFixed(4)}} Da</div>
                    <div><strong>Target M/Z:</strong> ${{c.mz.toFixed(4)}}</div>
                    <div style="margin-top:5px; color:#64748b; font-size:11px; word-break:break-all;"><strong>SMILES:</strong> ${{c.smiles}}</div>
                `;
                let imgDiv = document.createElement('div'); imgDiv.className = 'struct-img';
                if(c.img) imgDiv.innerHTML = c.img;
                card.appendChild(info); card.appendChild(imgDiv); listContainer.appendChild(card);
            }});
        }}
        if(select.value) renderPlate();
    </script>
</body>
</html>
"""
    return html_template

# ==========================================
# 3. Main Operational Execution
# ==========================================
if uploaded_file is not None:
    with st.spinner("Processing structural entries and executing stride partitioning..."):
        
        local_tmp_path = "temp_upload_data.sdf"
        with open(local_tmp_path, "wb") as f:
            f.write(uploaded_file.getvalue())
            
        raw_df = process_sdf(local_tmp_path)
        if os.path.exists(local_tmp_path): os.remove(local_tmp_path)
        
        if not raw_df.empty:
            # Assign pools using custom variables
            final_map = assign_wells_advanced(
                raw_df, 
                target_size=pool_size, 
                prefix=plate_prefix, 
                vol_comp=vol_per_comp,
                assay_vol=dest_well_vol,
                assay_conc=desired_conc
            )
            
            # Sort securely
            final_map = final_map.sort_values(by=['Destination_Plate', 'Destination_Well', 'Well_Sub_Index']).reset_index(drop=True)
            html_raw_data = final_map.copy()
            
            # Analytical validation checks
            total_wells_created = len(final_map['Destination_Well'].unique())
            wells_with_remainder = final_map[final_map['Compounds_In_Pool'] < pool_size]['Destination_Well'].nunique()
            complete_wells = total_wells_created - wells_with_remainder
            
            violating_wells_df = final_map[final_map['Min_m_z_Delta_In_Well'] < min_mz_threshold]
            num_violations = violating_wells_df['Destination_Well'].nunique()
            
            # System Metrics Dashboard
            st.success("Library compilation successful!")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Total Compounds Ingested", len(raw_df))
            col2.metric("Complete Target Pools", complete_wells)
            col3.metric("Partial (Back-Flush) Pools", wells_with_remainder)
            col4.metric("Calculated Echo Transfer (nL)", f"{final_map['Echo_Transfer_To_Destination_nL'].iloc[0]} nL")
            
            # Display Mass Spectrometer Alerts
            if num_violations > 0:
                st.warning(f"⚠️ Mass Resolution Alert: {num_violations} well pools contain compounds with a mass separation delta below your preferred {min_mz_threshold} Da threshold. Review the data preview table to verify safety.")
            else:
                st.info(f"✅ Mass Resolution Verified: All well pools maintain an internal mass delta greater than {min_mz_threshold} Da.")
                
            # ==========================================
            # 4. Advanced Excel Workbook Export Construction
            # ==========================================
            st.markdown("### Download Plate Layouts & Deliverables")
            down_col1, down_col2 = st.columns(2)
            
            # Building native openpyxl Excel buffer
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                final_map.to_excel(writer, sheet_name="ASMS Pooling Manifest", index=False)
            excel_buffer.seek(0)
            
            down_col1.download_button(
                label="Download Compound Management Manifest (Excel Workbook)",
                data=excel_buffer,
                file_name=f"{plate_prefix.strip().lower()}_pooling_manifest.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True
            )
            
            html_content = generate_interactive_html(html_raw_data)
            down_col2.download_button(
                label="Download Interactive Visual Layout (HTML Map)",
                data=html_content,
                file_name=f"{plate_prefix.strip().lower()}_interactive_layout.html",
                mime="text/html",
                use_container_width=True
            )
            
            # Plate Preview Layout View
            st.markdown("### Processed Plate Layout Preview")
            st.dataframe(
                final_map[[
                    'Destination_Plate', 'Destination_Well', 'Well_Sub_Index', 'Compounds_In_Pool', 'Backflush_Required',
                    'NCGC_ID', 'Exact_Mass', 'Target_m_z', 'Min_m_z_Delta_In_Well', 
                    'DMSO_Backflush_Volume_nL', 'Echo_Transfer_To_Destination_nL', 'Ionization_Mode'
                ]], use_container_width=True
            )
        else:
            st.error("No chemical structures successfully parsed. Verify the structural field blocks inside your uploaded SDF file.")
