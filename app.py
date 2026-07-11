import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import os

st.set_page_config(page_title="AS-MS Pooling Engine", page_icon="🧪", layout="wide")

st.title("NCATS ASMS 384-Well Pooling Engine")
st.markdown("Created by Colin Gordy for the creation of pooled compound libraries for HRMS analysis. This tool supports the development of an NCATS HTS screening assay utilizing His-tagged proteins bound to Ni-NTA magnetic beads. Upload Spotfire `.sdf` library export with the NCGC ID and SMILES to predict ionization modes, maximize $m/z$ diversity, and generate 384-well plate maps for CoMa source plate prep.")

# 1. File Uploader & Custom Settings
uploaded_file = st.file_uploader("Choose an SDF file", type=["sdf"])
plate_prefix = st.text_input("Plate Name Prefix", value="ASMS", help="Customize the starting label for your destination plates (e.g., NPC_ASMS, MIPE, etc.)")

def process_sdf(file_path):
    supplier = Chem.SDMolSupplier(file_path)
    data = []
    
    # Structural SMARTS patterns for ionization
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
            # Strip batch numbers (e.g., NCGC00091454-07 -> NCGC00091454)
            if '-' in str(sample_id):
                sample_id = str(sample_id).split('-')[0]
            
        try:
            exact_mass = Descriptors.ExactMolWt(mol)
            smiles = Chem.MolToSmiles(mol)
            
            # Skip empty records with no structural atoms, mass, or SMILES
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
                predicted_mode = "positive"  # Small molecule default
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

def assign_384_wells(df, pool_size=10, prefix="ASMS"):
    pooled_records = []
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    columns = range(1, 25)
    well_coordinates = [f"{r}{c:02d}" for r in rows for c in columns]
    
    current_plate = 1
    well_pointer = 0
    
    clean_prefix = prefix.strip().rstrip('_')
    
    for mode, group in df.groupby('Predicted_Mode'):
        sorted_group = group.sort_values(by='Target_m_z').reset_index(drop=True)
        num_compounds = len(sorted_group)
        num_pools = num_compounds // pool_size if num_compounds // pool_size > 0 else 1
        
        stride_pools = [[] for _ in range(num_pools)]
        for idx, row in sorted_group.iterrows():
            pool_idx = idx % num_pools
            if len(stride_pools[pool_idx]) < pool_size:
                stride_pools[pool_idx].append(row)
        
        for pool in stride_pools:
            if not pool: continue
            if well_pointer >= len(well_coordinates):
                well_pointer = 0
                current_plate += 1
                
            assigned_well = well_coordinates[well_pointer]
            well_pointer += 1
            
            for sub_idx, comp in enumerate(pool):
                pooled_records.append({
                    'Destination_Plate': f"{clean_prefix}_PLT_{current_plate}",  # 🛠️ FIXED: Standardized unified plate tracking label
                    'Destination_Well': assigned_well,
                    'Well_Sub_Index': sub_idx + 1,
                    'SAMPLE_ID': comp['SAMPLE_ID'],
                    'Exact_Mass': round(comp['Exact_Mass'], 4),
                    'Target_m_z': round(comp['Target_m_z'], 4),
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
        
        if plt not in plate_data_dict:
            plate_data_dict[plt] = {}
        if well not in plate_data_dict[plt]:
            plate_data_dict[plt][well] = []
            
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
            'id': row['SAMPLE_ID'],
            'mass': row['Exact_Mass'],
            'mz': row['Target_m_z'],
            'smiles': row['SMILES'],
            'img': svg_text
        })
        
    js_data_payload = json.dumps(plate_data_dict)
    
    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>NCATS ASMS Interactive Plate Mapper</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; background-color: #f8f9fa; margin: 0; padding: 20px; color: #333; }}
        .header {{ display: flex; justify-space-between; align-items: center; border-bottom: 2px solid #e9ecef; padding-bottom: 15px; margin-bottom: 20px; }}
        h1 {{ margin: 0; font-size: 24px; color: #1e293b; }}
        select {{ padding: 8px 16px; font-size: 16px; border-radius: 6px; border: 1px solid #cbd5e1; outline: none; background: white; cursor: pointer; }}
        .main-container {{ display: flex; gap: 24px; align-items: flex-start; }}
        .plate-box {{ background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); border: 1px solid #e2e8f0; }}
        .grid-384 {{ display: grid; grid-template-columns: 30px repeat(24, 26px); gap: 4px; align-items: center; justify-items: center; }}
        .col-header {{ font-size: 11px; font-weight: bold; color: #64748b; text-align: center; width: 26px; }}
        .row-header {{ font-size: 11px; font-weight: bold; color: #64748b; text-align: center; height: 26px; display: flex; align-items: center; justify-content: center; }}
        .well {{ width: 22px; height: 22px; border-radius: 50%; background-color: #f1f5f9; border: 1px solid #cbd5e1; cursor: pointer; transition: all 0.15s ease; position: relative; }}
        .well:hover {{ transform: scale(1.2); border-color: #475569; box-shadow: 0 0 4px rgba(0,0,0,0.2); }}
        .well.active-well {{ border-color: #2563eb; background-color: #dbeafe !important; box-shadow: 0 0 0 2px #2563eb; }}
        .well.populated {{ background-color: #94a3b8; }}
        .display-panel {{ flex: 1; min-width: 400px; background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); border: 1px solid #e2e8f0; max-height: 85vh; overflow-y: auto; }}
        .panel-title {{ font-size: 18px; font-weight: bold; margin-bottom: 15px; color: #1e293b; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; }}
        .compound-card {{ display: flex; align-items: center; gap: 15px; padding: 12px; margin-bottom: 12px; border: 1px solid #e2e8f0; border-radius: 8px; background: #f8fafc; }}
        .compound-info {{ flex: 1; font-size: 13px; }}
        .compound-id {{ font-size: 15px; font-weight: bold; color: #2563eb; margin-bottom: 4px; }}
        .struct-img {{ width: 130px; height: 130px; background: white; border: 1px solid #e2e8f0; border-radius: 6px; display: flex; align-items: center; justify-content: center; overflow: hidden; }}
        .struct-img img, .struct-img svg {{ max-width: 100%; max-height: 100%; }}
        .placeholder-text {{ color: #94a3b8; font-style: italic; text-align: center; margin-top: 50px; }}
    </style>
</head>
<body>

    <div class="header" style="display: flex; justify-content: space-between; width: 100%;">
        <h1>NCATS ASMS Interactive Pool Navigator</h1>
        <select id="plateSelect" onchange="renderPlate()"></select>
    </div>

    <div class="main-container">
        <div class="plate-box">
            <div id="gridContainer" class="grid-384"></div>
        </div>
        <div class="display-panel">
            <div id="panelTitle" class="panel-title">Well Pool Inspector</div>
            <div id="compoundsContainer">
                <div class="placeholder-text">Click any populated well on the left layout grid to inspect its 10 pooled compound structures...</div>
            </div>
        </div>
    </div>

    <script>
        const db = {js_data_payload};
        const rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'];
        
        const select = document.getElementById('plateSelect');
        Object.keys(db).sort().forEach(plt => {{
            let opt = document.createElement('option');
            opt.value = plt;
            opt.innerHTML = plt;
            select.appendChild(opt);
        }});

        function renderPlate() {{
            const currentPlate = select.value;
            const container = document.getElementById('gridContainer');
            container.innerHTML = '';
            
            let spacer = document.createElement('div');
            container.appendChild(spacer);
            
            for(let c=1; c<=24; c++) {{
                let header = document.createElement('div');
                header.className = 'col-header';
                header.innerHTML = c < 10 ? '0'+c : c;
                container.appendChild(header);
            }}
            
            rows.forEach(r => {{
                let rHeader = document.createElement('div');
                rHeader.className = 'row-header';
                rHeader.innerHTML = r;
                container.appendChild(rHeader);
                
                for(let c=1; c<=24; c++) {{
                    let wellName = r + (c < 10 ? '0'+c : c);
                    let wellDiv = document.createElement('div');
                    wellDiv.className = 'well';
                    wellDiv.id = wellName;
                    
                    const dynamicData = db[currentPlate] && db[currentPlate][wellName];
                    if (dynamicData) {{
                        wellDiv.classList.add('populated');
                        wellDiv.title = wellName + " (" + dynamicData.length + " compounds)";
                        wellDiv.onclick = () => selectWell(wellName, dynamicData, wellDiv);
                    }} else {{
                        wellDiv.title = wellName + " (Empty)";
                    }}
                    container.appendChild(wellDiv);
                }}
            }});
        }}

        function selectWell(wellName, compounds, element) {{
            document.querySelectorAll('.well').forEach(w => w.classList.remove('active-well'));
            element.classList.add('active-well');
            
            document.getElementById('panelTitle').innerHTML = "Contents of Well: " + wellName + " (" + compounds.length + " Compounds)";
            const listContainer = document.getElementById('compoundsContainer');
            listContainer.innerHTML = '';
            
            compounds.forEach(c => {{
                let card = document.createElement('div');
                card.className = 'compound-card';
                
                let info = document.createElement('div');
                info.className = 'compound-info';
                info.innerHTML = `
                    <div class="compound-id">${{c.id}}</div>
                    <div><strong>Exact Mass:</strong> ${{c.mass.toFixed(4)}} Da</div>
                    <div><strong>Target M/Z:</strong> ${{c.mz.toFixed(4)}}</div>
                    <div style="margin-top:5px; color:#64748b; font-size:11px; word-break:break-all;"><strong>SMILES:</strong> ${{c.smiles}}</div>
                `;
                
                let imgDiv = document.createElement('div');
                imgDiv.className = 'struct-img';
                if(c.img && c.img.trim() !== "") {{
                    imgDiv.innerHTML = c.img;
                }} else {{
                    imgDiv.innerHTML = `<span style="color:#cbd5e1; font-size:11px;">No Structure</span>`;
                }}
                
                card.appendChild(info);
                card.appendChild(imgDiv);
                listContainer.appendChild(card);
            }});
        }}

        if(select.value) renderPlate();
    </script>
</body>
</html>
"""
    return html_template

# ==========================================
# 2. Main App Execution
# ==========================================
if uploaded_file is not None:
    with st.spinner("Parsing chemical structures and generating optimal layout..."):
        
        local_tmp_path = "temp_upload_data.sdf"
        with open(local_tmp_path, "wb") as f:
            f.write(uploaded_file.getvalue())
            
        raw_df = process_sdf(local_tmp_path)
        
        if os.path.exists(local_tmp_path):
            os.remove(local_tmp_path)
        
        if not raw_df.empty:
            final_map = assign_384_wells(raw_df, pool_size=10, prefix=plate_prefix)
            final_map = final_map.sort_values(by=['Destination_Plate', 'Target_m_z']).reset_index(drop=True)
            
            html_raw_data = final_map.copy()
            
            final_map = final_map.rename(columns={
                'SAMPLE_ID': 'NCGC_ID',
                'Target_m_z': 'Target_M/Z',
                'Exact_Mass': 'Exact Mass'
            })
            
            st.success(f"Successfully processed {len(raw_df)} compounds!")
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Compounds", len(raw_df))
            col2.metric("Total 384-Well Pools Created", len(final_map['Destination_Well'].unique()))
            col3.metric("Total Plates Required", len(final_map['Destination_Plate'].unique()))
            
            # 3. Download Section
            st.markdown("### Download Plate Maps & Layouts")
            down_col1, down_col2 = st.columns(2)
            
            csv_data = final_map.to_csv(index=False).encode('utf-8')
            down_col1.download_button(
                label="Download Compound Management Manifest (CSV)",
                data=csv_data,
                file_name=f"{plate_prefix.strip().lower()}_pooling_manifest.csv",
                mime="text/csv",
                use_container_width=True
            )
            
            html_content = generate_interactive_html(html_raw_data)
            down_col2.download_button(
                label="Download Interactive HTML Plate Map (Visual Layout)",
                data=html_content,
                file_name=f"{plate_prefix.strip().lower()}_interactive_layout.html",
                mime="text/html",
                use_container_width=True
            )
            
            # 4. Interactive Preview
            st.markdown("### Plate Layout Preview")
            st.dataframe(
                final_map[['Destination_Plate', 'Destination_Well', 'Well_Sub_Index', 'NCGC_ID', 'Exact Mass', 'Target_M/Z', 'Ionization_Mode']], 
                use_container_width=True
            )
        else:
            st.error("Could not parse any valid structures from the uploaded SDF file.")
