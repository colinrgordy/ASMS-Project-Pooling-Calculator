import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import os
import io
import math
import json

st.set_page_config(page_title="ASMS Compound Suite", page_icon="⚜", layout="wide")

st.title("NCATS ASMS Compound Pooling & Quality Control Suite")

# Top Navigation Tabs
tab1, tab2, tab3, tab4 = st.tabs([
    "⚜ 1. Pooling Engine", 
    "🧪 2. Map Unpivoter", 
    "📊 3. Survey Pre-Filter", 
    "⚡ 4. Post-Run Reconciler"
])

# ==========================================
# TAB 1: MAIN POOLING ENGINE
# ==========================================
with tab1:
    st.markdown(
        "Created by Colin Gordy for use in the development of a semi-automated, small-molecule "
        "binders discovery assay utilizing HRMS. Compiles library entries into consolidated 384-well acoustic "
        "source pools, tracks volume normalization, and maps nanoliter transfers to 96-well target plates."
    )

    st.sidebar.header("Configuration Panel")

    st.sidebar.subheader("1. Library Pooling Options")
    pool_size = st.sidebar.number_input("Target Compounds per Well", min_value=2, max_value=50, value=10, step=1)
    min_mz_threshold = st.sidebar.number_input("Minimum Allowed Δm/z Threshold (Da)", min_value=0.5, max_value=10.0, value=2.0, step=0.5)

    st.sidebar.subheader("2. Volumetric Normalization & Physics")
    lib_stock_conc = st.sidebar.number_input("1536 Library Stock Concentration (mM)", min_value=1.0, max_value=100.0, value=10.0, step=1.0)
    vol_per_comp = st.sidebar.number_input("Source Plate: Aliquot Vol per Compound (nL)", min_value=10, max_value=5000, value=1000, step=100, help="Volume of each compound pipetted from the 1536 plate into the 384 source well.")
    target_source_vol_ul = st.sidebar.number_input("Source Plate: Target Total Well Volume (µL)", min_value=2.0, max_value=50.0, value=10.0, step=1.0, help="The desired final working volume inside the 384 Echo source well (typically 8-10 µL).")

    st.sidebar.subheader("3. Echo Assay Calculator")
    dest_well_vol = st.sidebar.number_input("Assay Plate: Total Well Volume (µL)", min_value=1.0, max_value=500.0, value=50.0, step=5.0)
    desired_conc = st.sidebar.number_input("Assay Plate: Target Concentration (µM)", min_value=0.1, max_value=100.0, value=10.0, step=1.0)

    plate_prefix = st.text_input(
        "Required: Plate Name Prefix", 
        value="", 
        placeholder="Ex: ASMS_NPC, ASMS_MGL, ASMS_MIPE, etc.",
        help="A unique identifier is strictly required to label data tracking tables and prevent accidental workbook file naming overlaps."
    )

    st.markdown("### Upload Core Campaign Assets")
    up_col1, up_col2 = st.columns(2)

    with up_col1:
        uploaded_file = st.file_uploader("Required: Choose an SDF Library File", type=["sdf"])

    with up_col2:
        st.info("💡 **Have a 2D visual map or volume survey?** Use Tabs 2 or 3 above to unpivot or pre-filter depleted wells first!")
        uploaded_inventory = st.file_uploader(
            "Optional: Upload 1536 Master Plate Maps", 
            type=["csv", "xlsx"], 
            help="Provide the manifest file containing real-world freezer locations to generate the initial 1536 to 384 pool picklist file."
        )

    def process_sdf(file_path):
        supplier = Chem.SDMolSupplier(file_path)
        data = []
        omissions = []
        
        basic_nitrogen = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(N-[#6a])]")
        acidic_group = Chem.MolFromSmarts("[C,S](=[O,S])[O;H1,-1]")
        
        for idx, mol in enumerate(supplier):
            if mol is None:
                omissions.append({
                    'SDF_Record_Index': idx + 1,
                    'NCGC_ID': "UNPARSEABLE_RECORD",
                    'Omission_Reason': "Corrupted / Dead SDF Record Block"
                })
                continue
            
            sample_id = None
            for prop_name in ['SAMPLE_ID', 'Name', 'ID', 'sample_id', 'id']:
                if mol.HasProp(prop_name):
                    sample_id = mol.GetProp(prop_name)
                    break
            
            if not sample_id: 
                sample_id = f"UNKNOWN_ID_REC_{idx + 1}"
            else:
                if '-' in str(sample_id):
                    sample_id = str(sample_id).split('-')[0]
            
            working_mol = mol
            if len(Chem.GetMolFrags(mol)) > 1:
                try:
                    frags = Chem.GetMolFrags(mol, asMols=True)
                    working_mol = max(frags, key=lambda m: m.GetNumAtoms())
                except:
                    pass
                
            try:
                exact_mass = Descriptors.ExactMolWt(working_mol)
                smiles = Chem.MolToSmiles(working_mol)
                
                if exact_mass == 0:
                    omissions.append({
                        'SDF_Record_Index': idx + 1,
                        'NCGC_ID': sample_id,
                        'Omission_Reason': "Empty Record (Zero Molecular Weight)"
                    })
                    continue
                
                if not smiles or smiles.strip() == "":
                    omissions.append({
                        'SDF_Record_Index': idx + 1,
                        'NCGC_ID': sample_id,
                        'Omission_Reason': "Missing Structural SMILES Data String"
                    })
                    continue
                    
                has_base = working_mol.HasSubstructMatch(basic_nitrogen)
                has_acid = working_mol.HasSubstructMatch(acidic_group)
                
                if has_base and not has_acid:
                    predicted_mode = "positive"
                    target_mz = exact_mass + 1.0073
                elif has_acid and not has_base:
                    predicted_mode = "negative"
                    target_mz = exact_mass - 1.0073
                else:
                    predicted_mode = "positive"  
                    target_mz = exact_mass + 1.0073
                    
                data.append({
                    'SAMPLE_ID': sample_id,
                    'Exact_Mass': exact_mass,
                    'Target_m_z': target_mz,
                    'Predicted_Mode': predicted_mode,
                    'SMILES': smiles
                })
            except Exception as e:
                omissions.append({
                    'SDF_Record_Index': idx + 1,
                    'NCGC_ID': sample_id,
                    'Omission_Reason': f"Unexpected Processing Exception: {str(e)}"
                })
                continue
                
        return pd.DataFrame(data), pd.DataFrame(omissions)

    def assign_wells_advanced(df, target_size, prefix, vol_comp, target_source_vol_ul, lib_stock_uM, assay_vol, assay_conc):
        pooled_records = []
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
        columns = range(1, 25) 
        well_coordinates = [f"{r}{c:02d}" for r in rows for c in columns]
        
        current_plate = 1
        well_pointer = 0
        clean_prefix = prefix.strip().rstrip('_')
        
        target_source_vol_nl = target_source_vol_ul * 1000.0
        source_well_conc_uM = lib_stock_uM * (vol_comp / target_source_vol_nl)
        
        assay_vol_nl = assay_vol * 1000.0
        echo_transfer_volume_nl = (assay_conc * assay_vol_nl) / source_well_conc_uM
        
        for mode, group in df.groupby('Predicted_Mode'):
            sorted_group = group.sort_values(by='Target_m_z').reset_index(drop=True)
            num_compounds = len(sorted_group)
            
            if num_compounds == 0: continue
            
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
                total_compound_fluid_nl = actual_pool_count * vol_comp
                dmso_backflush_nl = target_source_vol_nl - total_compound_fluid_nl
                
                pool_mzs = sorted([comp['Target_m_z'] for comp in pool])
                if len(pool_mzs) > 1:
                    min_delta_observed = min(pool_mzs[i+1] - pool_mzs[i] for i in range(len(pool_mzs)-1))
                else:
                    min_delta_observed = float('inf')
                    
                for sub_idx, comp in enumerate(pool):
                    pooled_records.append({
                        'Source_Plate_384': f"{clean_prefix}_SRC_PLT_{current_plate}",
                        'Source_Well_384': assigned_well,
                        'Well_Sub_Index': sub_idx + 1,
                        'Compounds_In_Pool': actual_pool_count,
                        'Backflush_Required': "YES" if dmso_backflush_nl > 0 else "NO",
                        'NCGC_ID': comp['SAMPLE_ID'],
                        'Exact_Mass': round(comp['Exact_Mass'], 4),
                        'Target_m_z': round(comp['Target_m_z'], 4),
                        'Min_Δm/z_In_Well': round(min_delta_observed, 4) if min_delta_observed != float('inf') else 0.0,
                        'DMSO_Backflush_Volume_nL': max(0.0, dmso_backflush_nl),
                        'Total_Well_Fluid_Vol_nL': target_source_vol_nl,
                        'Assay_Total_Volume_µL': assay_vol,
                        'Assay_Target_Conc_µM': assay_conc,
                        'Echo_Transfer_Volume_nL': round(echo_transfer_volume_nl, 2),
                        'Ionization_Mode': comp['Predicted_Mode'],
                        'SMILES': comp['SMILES']
                    })
                    
        return pd.DataFrame(pooled_records)

    def generate_dual_interactive_html(df, target_pool_max):
        source_dict = {}
        assay_dict = {}
        
        for _, row in df.iterrows():
            src_plt = row['Source_Plate_384']
            src_well = row['Source_Well_384']
            asy_plt = row['Assay_Plate_96']
            asy_well = row['Assay_Well_96']
            
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
                
            comp_card = {
                'id': row['NCGC_ID'],
                'mass': row['Exact_Mass'],
                'mz': row['Target_m_z'],
                'smiles': row['SMILES'],
                'img': svg_text,
                'mode': row['Ionization_Mode'],
                'backflush': int(row['DMSO_Backflush_Volume_nL']),
                'actual_count': int(row['Compounds_In_Pool']),
                'target_count': int(target_pool_max)
            }
            
            if src_plt not in source_dict: source_dict[src_plt] = {}
            if src_well not in source_dict[src_plt]: source_dict[src_plt][src_well] = []
            source_dict[src_plt][src_well].append(comp_card)
            
            if asy_plt not in assay_dict: assay_dict[asy_plt] = {}
            if asy_well not in assay_dict[asy_plt]: assay_dict[asy_plt][asy_well] = []
            assay_dict[asy_plt][asy_well].append(comp_card)
            
        full_payload = json.dumps({'source': source_dict, 'assay': assay_dict})
        
        html_template = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>NCATS Dual-Plate Assay Navigator</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; background-color: #f8f9fa; margin: 0; padding: 20px; color: #333; }
        .header { display: flex; justify-content: space-between; align-items: center; border-bottom: 2px solid #e9ecef; padding-bottom: 15px; margin-bottom: 20px; gap: 15px; flex-wrap: wrap; }
        .controls { display: flex; gap: 12px; align-items: center; }
        h1 { margin: 0; font-size: 22px; color: #1e293b; }
        select, button { padding: 8px 16px; font-size: 14px; border-radius: 6px; border: 1px solid #cbd5e1; background: white; cursor: pointer; font-weight: 500; }
        button.active-btn { background-color: #2563eb; color: white; border-color: #2563eb; }
        
        .wrapper-nav { display: flex; align-items: center; gap: 0px; background: white; border: 1px solid #cbd5e1; border-radius: 6px; padding: 0px; overflow: hidden; }
        .wrapper-nav select { border: none; border-radius: 0; padding: 8px 12px; background: transparent; outline: none; }
        .btn-arrow { border: none; background: transparent; padding: 8px 12px; font-size: 11px; color: #64748b; transition: all 0.1s ease; border-radius: 0; }
        .btn-arrow:hover { background-color: #f1f5f9; color: #1e293b; }
        .btn-arrow:first-child { border-right: 1px solid #e2e8f0; }
        .btn-arrow:last-child { border-left: 1px solid #e2e8f0; }
        
        .main-container { display: flex; gap: 24px; align-items: flex-start; }
        .plate-box { background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); border: 1px solid #e2e8f0; }
        .map-legend { display: flex; gap: 24px; margin-bottom: 20px; font-size: 13px; font-weight: 600; justify-content: center; background: #f8fafc; padding: 10px; border-radius: 8px; border: 1px solid #e2e8f0; flex-wrap: wrap; }
        .legend-item { display: flex; align-items: center; gap: 8px; }
        
        .grid-container { display: grid; gap: 4px; align-items: center; justify-items: center; }
        .grid-384 { grid-template-columns: 30px repeat(24, 26px); }
        .grid-96 { grid-template-columns: 30px repeat(12, 40px); }
        
        .col-header { font-size: 11px; font-weight: bold; color: #64748b; text-align: center; }
        .row-header { font-size: 11px; font-weight: bold; color: #64748b; text-align: center; display: flex; align-items: center; justify-content: center; }
        
        .well { border-radius: 50%; background-color: #f1f5f9; border: 1px solid #cbd5e1; cursor: pointer; transition: all 0.15s ease; box-sizing: border-box; }
        .well-384 { width: 22px; height: 22px; }
        .well-96 { width: 34px; height: 34px; }
        
        .well.populated.positive { background-color: #bfdbfe; border-color: #3b82f6; }
        .well.populated.negative { background-color: #fecdd3; border-color: #f43f5e; }
        .well.backflush-needed { border-style: dashed !important; border-width: 2px !important; }
        .well.incomplete-pool { border-style: dashed !important; border-width: 2px !important; border-color: #ea580c !important; }
        .well:hover { transform: scale(1.15); border-color: #475569 !important; box-shadow: 0 0 4px rgba(0,0,0,0.15); z-index: 10; }
        .well.active-well { border-color: #1e3a8a !important; background-color: #eff6ff !important; box-shadow: 0 0 0 3px #3b82f6; }
        
        .display-panel { flex: 1; min-width: 400px; background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); border: 1px solid #e2e8f0; max-height: 85vh; overflow-y: auto; }
        .panel-title { font-size: 18px; font-weight: bold; margin-bottom: 15px; color: #1e293b; border-bottom: 1px solid #e2e8f0; padding-bottom: 8px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 10px; }
        .backflush-tag { font-size: 12px; font-weight: bold; color: #b45309; background-color: #fef3c7; padding: 4px 10px; border-radius: 6px; border: 1px solid #fde68a; }
        .warning-tag { font-size: 12px; font-weight: bold; color: #dc2626; background-color: #fef2f2; padding: 4px 10px; border-radius: 6px; border: 1px solid #fca5a5; }
        
        .compound-card { display: flex; align-items: center; gap: 15px; padding: 12px; margin-bottom: 12px; border: 1px solid #e2e8f0; border-radius: 8px; background: #f8fafc; }
        .compound-info { flex: 1; font-size: 13px; }
        .compound-id { font-size: 15px; font-weight: bold; color: #2563eb; margin-bottom: 4px; }
        .struct-img { width: 130px; height: 130px; background: white; border: 1px solid #e2e8f0; border-radius: 6px; display: flex; align-items: center; justify-content: center; overflow: hidden; }
        .struct-img svg { max-width: 100%; max-height: 100%; }
        .placeholder-text { color: #94a3b8; font-style: italic; text-align: center; margin-top: 50px; }
    </style>
</head>
<body>

    <div class="header">
        <h1>ASMS Project Interactive Navigator</h1>
        <div class="controls">
            <button id="btn384" class="active-btn" onclick="setViewType('source')">384-Well Source Plates</button>
            <button id="btn96" onclick="setViewType('assay')">96-Well Assay Plates</button>
            
            <div class="wrapper-nav">
                <button class="btn-arrow" onclick="pagePlates(-1)" title="Previous Plate">◀</button>
                <select id="plateSelect" onchange="renderPlate()"></select>
                <button class="btn-arrow" onclick="pagePlates(1)" title="Next Plate">▶</button>
            </div>
        </div>
    </div>

    <div class="main-container">
        <div class="plate-box">
            <div id="legendBox" class="map-legend"></div>
            <div id="gridContainer" class="grid-container"></div>
        </div>
        <div class="display-panel">
            <div id="panelTitle" class="panel-title">Well Pool Inspector</div>
            <div id="compoundsContainer">
                <div class="placeholder-text">Click any populated well on the left layout grid to inspect its contents...</div>
            </div>
        </div>
    </div>

    <script>
        const masterDataset = {js_data_payload};
        let currentViewMode = 'source'; 
        
        const rows384 = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'];
        const rows96 = ['A','B','C','D','E','F','G','H'];
        
        function setViewType(mode) {
            currentViewMode = mode;
            document.getElementById('btn384').classList.toggle('active-btn', mode === 'source');
            document.getElementById('btn96').classList.toggle('active-btn', mode === 'assay');
            
            const select = document.getElementById('plateSelect');
            select.innerHTML = '';
            const targetedSubTree = masterDataset[currentViewMode];
            
            Object.keys(targetedSubTree).sort().forEach(plt => {
                let opt = document.createElement('option');
                opt.value = plt; opt.innerHTML = plt; select.appendChild(opt);
            });
            
            renderLegend();
            renderPlate();
        }

        function pagePlates(direction) {
            const select = document.getElementById('plateSelect');
            const proposedIndex = select.selectedIndex + direction;
            if (proposedIndex >= 0 && proposedIndex < select.options.length) {
                select.selectedIndex = proposedIndex;
                renderPlate();
            }
        }

        function renderLegend() {
            const legend = document.getElementById('legendBox');
            if (currentViewMode === 'source') {
                legend.innerHTML = `
                    <div class="legend-item"><div style="width:14px; height:14px; border-radius:50%; background-color:#fecdd3; border:1px solid #f43f5e;"></div><span>Negative Pools</span></div>
                    <div class="legend-item"><div style="width:14px; height:14px; border-radius:50%; background-color:#bfdbfe; border:1px solid #3b82f6;"></div><span>Positive Pools</span></div>
                    <div class="legend-item"><div style="width:14px; height:14px; border-radius:50%; background-color:#f1f5f9; border:2px dashed #475569;"></div><span>Dashed Border = Requires DMSO Back-flush</span></div>
                `;
            } else {
                legend.innerHTML = `
                    <div class="legend-item"><div style="width:14px; height:14px; border-radius:50%; background-color:#fecdd3; border:1px solid #f43f5e;"></div><span>Negative Assay Well Block</span></div>
                    <div class="legend-item"><div style="width:14px; height:14px; border-radius:50%; background-color:#bfdbfe; border:1px solid #3b82f6;"></div><span>Positive Assay Well Block</span></div>
                    <div class="legend-item"><div style="width:14px; height:14px; border-radius:50%; background-color:#f1f5f9; border:2px dashed #ea580c;"></div><span>Orange Dashed Border = Incomplete Pool (Fewer Peaks)</span></div>
                `;
            }
        }

        function renderPlate() {
            const select = document.getElementById('plateSelect');
            const currentPlate = select.value;
            const container = document.getElementById('gridContainer');
            container.innerHTML = '';
            
            if(!currentPlate) return;
            
            const is384 = currentViewMode === 'source';
            container.className = is384 ? "grid-container grid-384" : "grid-container grid-96";
            
            let spacer = document.createElement('div'); container.appendChild(spacer);
            
            const numCols = is384 ? 24 : 12;
            const targetRows = is384 ? rows384 : rows96;
            
            for(let c=1; c<=numCols; c++) {
                let header = document.createElement('div'); header.className = 'col-header';
                header.innerHTML = c < 10 ? '0'+c : c; container.appendChild(header);
            }
            
            targetRows.forEach(r => {
                let rHeader = document.createElement('div'); rHeader.className = 'row-header';
                rHeader.innerHTML = r; container.appendChild(rHeader);
                
                for(let c=1; c<=numCols; c++) {
                    let wellName = r + (c < 10 ? '0'+c : c);
                    let wellDiv = document.createElement('div');
                    wellDiv.className = is384 ? 'well well-384' : 'well well-96';
                    wellDiv.id = wellName;
                    
                    const dynamicData = masterDataset[currentViewMode][currentPlate] && masterDataset[currentViewMode][currentPlate][wellName];
                    if (dynamicData && dynamicData.length > 0) {
                        wellDiv.classList.add('populated');
                        const wellMode = dynamicData[0].mode;
                        wellDiv.classList.add(wellMode);
                        
                        if (is384 && dynamicData[0].backflush > 0) {
                            wellDiv.classList.add('backflush-needed');
                        }
                        if (!is384 && dynamicData[0].actual_count < dynamicData[0].target_count) {
                            wellDiv.classList.add('incomplete-pool');
                        }
                        wellDiv.onclick = () => selectWell(wellName, dynamicData, wellDiv);
                    }
                    container.appendChild(wellDiv);
                }
            });
        }

        function selectWell(wellName, compounds, element) {
            document.querySelectorAll('.well').forEach(w => w.classList.remove('active-well'));
            element.classList.add('active-well');
            
            const wellModeLabel = compounds[0].mode.toUpperCase();
            const plateContextLabel = currentViewMode === 'source' ? 'Source Well' : 'Assay Well';
            
            let headerHTML = `<span>Contents of ${plateContextLabel}: ${wellName} (${wellModeLabel} Mode)</span>`;
            if(currentViewMode === 'source' && compounds[0].backflush > 0) {
                headerHTML += `<span class="backflush-tag">⚠️ DMSO Back-flush: +${compounds[0].backflush} nL</span>`;
            }
            if(currentViewMode === 'assay' && compounds[0].actual_count < compounds[0].target_count) {
                headerHTML += `<span class="warning-tag">⚠️ Incomplete Pool: ${compounds[0].actual_count}/${compounds[0].target_count} Compounds</span>`;
            }
            document.getElementById('panelTitle').innerHTML = headerHTML;
            
            const listContainer = document.getElementById('compoundsContainer');
            listContainer.innerHTML = '';
            
            compounds.forEach(c => {
                let card = document.createElement('div'); card.className = 'compound-card';
                let info = document.className = 'compound-info';
                card.innerHTML = `
                    <div style="flex:1;">
                        <div style="font-size:15px; font-weight:bold; color:#2563eb; margin-bottom:4px;">${c.id}</div>
                        <div><strong>Exact Mass:</strong> ${c.mass.toFixed(4)} Da</div>
                        <div><strong>Target M/Z:</strong> ${c.mz.toFixed(4)}</div>
                        <div style="margin-top:5px; color:#64748b; font-size:11px; word-break:break-all;"><strong>SMILES:</strong> ${c.smiles}</div>
                    </div>
                    <div style="width:130px; height:130px; background:white; border:1px solid #e2e8f0; border-radius:6px; display:flex; align-items:center; justify-content:center; overflow:hidden;">
                        ${c.img ? c.img : ''}
                    </div>
                `;
                listContainer.appendChild(card);
            });
        }

        setViewType('source');
    </script>
</body>
</html>
"""
        return html_template.replace("{js_data_payload}", full_payload)

    # Main Pipeline Execution
    if uploaded_file is not None:
        if not plate_prefix.strip():
            st.error("⚠️ **Missing Required Field:** Enter a unique Plate Name Prefix above before running calculations.")
            st.stop()

        max_possible_compound_vol_nl = pool_size * vol_per_comp
        target_source_vol_nl = target_source_vol_ul * 1000.0
        
        if max_possible_compound_vol_nl > target_source_vol_nl:
            st.error(f"❌ **Physical Fluidic Paradox Error:** You requested a pool size of **{pool_size} compounds** at **{vol_per_comp} nL** each. This requires **{max_possible_compound_vol_nl / 1000.0} µL** per well from library aliquots alone, physically overflowing your target source well capacity of **{target_source_vol_ul} µL**.")
            st.stop()

        with st.spinner("Executing double-stage library calculations..."):
            local_tmp_path = "temp_upload_data.sdf"
            with open(local_tmp_path, "wb") as f:
                f.write(uploaded_file.getvalue())
                
            raw_df, skipped_df = process_sdf(file_path=local_tmp_path)
            if os.path.exists(local_tmp_path): os.remove(local_tmp_path)
            
            if not raw_df.empty:
                lib_stock_uM = lib_stock_conc * 1000.0
                
                source_map = assign_wells_advanced(
                    raw_df, 
                    target_size=pool_size, 
                    prefix=plate_prefix, 
                    vol_comp=vol_per_comp,
                    target_source_vol_ul=target_source_vol_ul,
                    lib_stock_uM=lib_stock_uM,
                    assay_vol=dest_well_vol,
                    assay_conc=desired_conc
                )
                
                unique_source_wells = source_map[['Source_Plate_384', 'Source_Well_384']].drop_duplicates().reset_index(drop=True)
                assay_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
                assay_cols = range(1, 13)
                assay_coordinates = [f"{r}{c:02d}" for r in assay_rows for c in assay_cols]
                
                coordinate_mapping_index = {}
                for idx, r_wells in unique_source_wells.iterrows():
                    plate_idx = (idx // 96) + 1
                    assigned_96_well = assay_coordinates[idx % 96]
                    coordinate_mapping_index[(r_wells['Source_Plate_384'], r_wells['Source_Well_384'])] = (f"{plate_prefix}_ASSAY_PLT_{plate_idx}", assigned_96_well)
                
                source_map['Assay_Plate_96'] = source_map.apply(lambda r: coordinate_mapping_index[(r['Source_Plate_384'], r['Source_Well_384'])][0], axis=1)
                source_map['Assay_Well_96'] = source_map.apply(lambda r: coordinate_mapping_index[(r['Source_Plate_384'], r['Source_Well_384'])][1], axis=1)
                source_map['Designated_Pool_Size'] = pool_size
                source_map['Actual_Pool_Size'] = source_map['Compounds_In_Pool']
                source_map['Pool_Status'] = source_map.apply(lambda r: "COMPLETE" if r['Actual_Pool_Size'] == pool_size else f"⚠️ INCOMPLETE ({r['Actual_Pool_Size']}/{pool_size})", axis=1)
                source_map = source_map.sort_values(by=['Source_Plate_384', 'Source_Well_384', 'Well_Sub_Index']).reset_index(drop=True)
                
                total_384_wells = len(source_map['Source_Well_384'].unique())
                total_96_wells = len(source_map[['Assay_Plate_96', 'Assay_Well_96']].drop_duplicates())
                plates_needed = math.ceil(total_96_wells / 96)
                backflush_wells_count = source_map[source_map['DMSO_Backflush_Volume_nL'] > 0]['Source_Well_384'].nunique()
                
                st.success("HTS Screening manifests successfully generated!")
                
                dash_col1, dash_col2, dash_col3, dash_col4 = st.columns(4)
                dash_col1.metric("Total Library Compounds", len(raw_df))
                dash_col2.metric("384-Well Source Pools", total_384_wells)
                dash_col3.metric("96-Well Assay Plates", plates_needed)
                dash_col4.metric("DMSO Back-Flush Actions", backflush_wells_count)
                
                with st.expander("Click to view structural omissions and data anomalies"):
                    if not skipped_df.empty:
                        st.warning(f"Total entries filtered out from raw input file: {len(skipped_df)}")
                        st.dataframe(skipped_df, use_container_width=True)
                    else:
                        st.info("Pristine Library File: 0 structural exclusions or parsing failures recorded.")
                
                violating_pools = source_map[source_map['Min_Δm/z_In_Well'] < min_mz_threshold]['Source_Well_384'].nunique()
                if violating_pools > 0:
                    st.warning(f"⚠️ Mass Resolution Alert: {violating_pools} well pools contain compounds falling below your preferred {min_mz_threshold} Da Δm/z resolution limit.")
                else:
                    st.info(f"✅ Mass Resolution Checked: All pooling mixtures maintain structural Δm/z separation limits above {min_mz_threshold} Da.")
                    
                # ==========================================
                # 5. Multi-Deliverable Export Hub
                # ==========================================
                st.markdown("### Download Campaign Assets")
                
                # Check for uploaded 1536 Master Map to generate Picklist #0
                if uploaded_inventory is not None:
                    try:
                        if uploaded_inventory.name.endswith('.csv'):
                            inv_df = pd.read_csv(uploaded_inventory)
                        else:
                            inv_df = pd.read_excel(uploaded_inventory)
                        
                        inv_df.columns = [str(c).strip() for c in inv_df.columns]
                        
                        expected_id_col = next((c for c in inv_df.columns if c.upper() in ['NCGC_ID', 'SAMPLE_ID', 'ID', 'COMPOUND_ID']), None)
                        expected_plt_col = next((c for c in inv_df.columns if 'PLATE' in c.upper() or 'BARCODE' in c.upper()), None)
                        expected_wel_col = next((c for c in inv_df.columns if 'WELL' in c.upper() or 'COORD' in c.upper()), None)
                        
                        if expected_id_col and expected_plt_col and expected_wel_col:
                            # 🛠️ FIXED: Normalize match_id on BOTH sides (strip batch suffix e.g., 'NCGC00015716-09' -> 'NCGC00015716')
                            inv_df['match_id'] = inv_df[expected_id_col].astype(str).str.strip().apply(lambda x: x.split('-')[0])
                            source_map['match_id'] = source_map['NCGC_ID'].astype(str).str.strip().apply(lambda x: x.split('-')[0])
                            
                            consolidation_df = pd.merge(source_map, inv_df, on='match_id', how='inner')
                            
                            if not consolidation_df.empty:
                                picklist_1536_to_384 = consolidation_df[[
                                    expected_plt_col, expected_wel_col, 'Source_Plate_384', 'Source_Well_384'
                                ]].copy()
                                picklist_1536_to_384['Transfer Volume'] = vol_per_comp
                                picklist_1536_to_384.columns = [
                                    'Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume'
                                ]
                                
                                picklist_1536_to_384 = picklist_1536_to_384.drop_duplicates().reset_index(drop=True)
                                
                                # 🛠️ SORT BY SOURCE PLATE FIRST TO ELIMINATE PLATE SWAPS
                                picklist_1536_to_384 = picklist_1536_to_384.sort_values(
                                    by=['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well']
                                ).reset_index(drop=True)
                                
                                buf_up = io.StringIO()
                                picklist_1536_to_384.to_csv(buf_up, index=False)
                                csv_up = buf_up.getvalue()
                                
                                st.success(f"✅ **1536 ➔ 384 Consolidation Picklist Generated!** Matched {len(consolidation_df)} compound locations.")
                                st.download_button(
                                    label="0. Download 1536 ➔ 384 Consolidation Picklist (.csv)",
                                    data=csv_up,
                                    file_name=f"{plate_prefix.strip().lower()}_1536_to_384_consolidation.csv",
                                    mime="text/csv",
                                    type="primary",
                                    use_container_width=False
                                )
                            else:
                                st.error("⚠️ **0 Match IDs Identified:** Make sure your 1536 Master Map contains an 'NCGC_ID' column matching the SDF.")
                        else:
                            st.error(f"⚠️ **Missing Required Mapping Columns:** Need 'NCGC_ID', 'Plate_1536', and 'Well_1536' in Master Map. Found columns: {list(inv_df.columns)}")
                    except Exception as ex:
                        st.error(f"Upstream pipeline error: {str(ex)}")

                down_col1, down_col2, down_col3, down_col4 = st.columns(4)

                src_excel_cols = [
                    'Source_Plate_384', 'Source_Well_384', 'Well_Sub_Index', 'Compounds_In_Pool', 'Backflush_Required',
                    'NCGC_ID', 'Exact_Mass', 'Target_m_z', 'Min_Δm/z_In_Well', 'DMSO_Backflush_Volume_nL', 'Total_Well_Fluid_Vol_nL'
                ]
                buf_src = io.BytesIO()
                with pd.ExcelWriter(buf_src, engine='openpyxl') as writer:
                    source_map[src_excel_cols].to_excel(writer, sheet_name="384_Source_Prep", index=False)
                buf_src.seek(0)
                down_col1.download_button(
                    label="1. Download Source Prep Workbook (.xlsx)",
                    data=buf_src,
                    file_name=f"{plate_prefix.strip().lower()}_source_prep_manifest.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
                
                asy_excel_cols = [
                    'Assay_Plate_96', 'Assay_Well_96', 'Pool_Status', 'Designated_Pool_Size', 'Actual_Pool_Size', 
                    'Source_Plate_384', 'Source_Well_384', 'NCGC_ID', 'Exact_Mass', 'Target_m_z', 
                    'Assay_Total_Volume_µL', 'Assay_Target_Conc_µM', 'Echo_Transfer_Volume_nL', 'Ionization_Mode'
                ]
                buf_asy = io.BytesIO()
                with pd.ExcelWriter(buf_asy, engine='openpyxl') as writer:
                    source_map[asy_excel_cols].to_excel(writer, sheet_name="96_Assay_Run", index=False)
                buf_asy.seek(0)
                down_col2.download_button(
                    label="2. Download Assay Run Manifest (.xlsx)",
                    data=buf_asy,
                    file_name=f"{plate_prefix.strip().lower()}_assay_run_manifest.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
                
                html_payload = generate_dual_interactive_html(df=source_map, target_pool_max=pool_size)
                down_col3.download_button(
                    label="3. Download Campaign Browser Map (.html)",
                    data=html_payload,
                    file_name=f"{plate_prefix.strip().lower()}_campaign_map.html",
                    mime="text/html",
                    use_container_width=True
                )
                
                echo_picklist_df = source_map[[
                    'Source_Plate_384', 'Source_Well_384', 'Assay_Plate_96', 'Assay_Well_96', 'Echo_Transfer_Volume_nL'
                ]].drop_duplicates().reset_index(drop=True)
                
                echo_picklist_df.columns = [
                    'Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume'
                ]
                
                buf_pick = io.StringIO()
                echo_picklist_df.to_csv(buf_pick, index=False)
                csv_pick = buf_pick.getvalue()
                
                down_col4.download_button(
                    label="4. Download Echo Run Picklist (.csv)",
                    data=csv_pick,
                    file_name=f"{plate_prefix.strip().lower()}_echo_assay_picklist.csv",
                    mime="text/csv",
                    use_container_width=True
                )
                
                st.markdown("### Unified Data Matrix Preview")
                st.dataframe(source_map[[
                    'Source_Plate_384', 'Source_Well_384', 'Well_Sub_Index', 'Backflush_Required',
                    'Assay_Plate_96', 'Assay_Well_96', 'Pool_Status', 'NCGC_ID', 'Exact_Mass', 'Target_m_z', 
                    'Min_Δm/z_In_Well', 'DMSO_Backflush_Volume_nL', 'Assay_Target_Conc_µM', 'Echo_Transfer_Volume_nL'
                ]], use_container_width=True)
                
            else:
                st.error("SDF parsing returned an empty data array. Check property tagging fields.")

# ==========================================
# TAB 2: PLATE MAP UNPIVOTER (CORRECTED)
# ==========================================
with tab2:
    st.subheader("🧪 Visual Plate Map Unpivoter")
    st.write("Convert 2D visual grid Excel sheets (A–AF rows) into flat CSV manifests.")

    uploaded_map_file = st.file_uploader("Upload Excel Plate Map (.xlsx)", type=["xlsx", "xls"], key="unpivoter_uploader")

    if uploaded_map_file is not None:
        try:
            xls_file = pd.ExcelFile(uploaded_map_file)

            all_records = []
            for sheet_name in xls_file.sheet_names:
                df_map = pd.read_excel(xls_file, sheet_name=sheet_name, header=None)
                for r_idx in range(df_map.shape[0]):
                    row_label = str(df_map.iloc[r_idx, 0]).strip()
                    if not row_label or row_label.lower() == 'nan': 
                        continue
                    
                    # 🛠️ FIXED: c_idx = 1 is 1536 Col 01 (Excel B), c_idx = 5 is 1536 Col 05 (Excel F)
                    for c_idx in range(1, df_map.shape[1]):
                        val = df_map.iloc[r_idx, c_idx]
                        if pd.notna(val) and str(val).strip().lower() != 'nan':
                            well_id = f"{row_label}{c_idx:02d}"  # Directly maps c_idx 5 -> 'A05'
                            all_records.append({
                                "NCGC_ID": str(val).strip(), 
                                "Plate_1536": sheet_name, 
                                "Well_1536": well_id
                            })
            
            flat_df = pd.DataFrame(all_records)
            
            st.metric("Total Linearized Wells Extracted", len(flat_df))
            st.dataframe(flat_df.head(20), use_container_width=True)
            
            buf = io.StringIO()
            flat_df.to_csv(buf, index=False)
            st.download_button("⬇️ Download Corrected Linearized CSV Map", buf.getvalue(), "1536_master_map_flat_corrected.csv", "text/csv", type="primary")
        except Exception as e:
            st.error(f"Error processing file: {e}")

# ==========================================
# TAB 3: ECHO SURVEY VOLUME PRE-FILTER
# ==========================================
with tab3:
    st.subheader("📊 Echo Survey Volume Pre-Filter")
    st.write("Cross-reference your 1536 master plate map against an Echo Volume Survey spreadsheet to filter out low-volume wells before running calculations.")

    col_s1, col_s2 = st.columns(2)
    with col_s1:
        s_map = st.file_uploader("1. Upload Linearized 1536 Master Map (.csv)", type=["csv"], key="s_map_up")
    with col_s2:
        s_survey = st.file_uploader("2. Upload Echo Survey Spreadsheet (.xlsx)", type=["xlsx", "xls"], key="s_surv_up")

    vol_cutoff = st.number_input("Minimum Dispense Volume Cutoff (µL)", min_value=0.5, max_value=5.0, value=1.5, step=0.1)

    if s_map is not None and s_survey is not None:
        try:
            map_df = pd.read_csv(s_map)
            df_survey_raw = pd.read_excel(s_survey, header=None)
            
            row_letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','AA','AB','AC','AD','AE','AF']
            
            blocks = []
            cur_b = []
            for r in range(len(df_survey_raw)):
                if len(df_survey_raw.iloc[r, 5:].dropna()) > 0: cur_b.append(r)
                else:
                    if cur_b: blocks.append(cur_b); cur_b = []
            if cur_b: blocks.append(cur_b)

            unique_map_plates = map_df['Plate_1536'].unique().tolist()
            survey_recs = []

            for p_idx, b in enumerate(blocks):
                plt_name = unique_map_plates[p_idx] if p_idx < len(unique_map_plates) else f"Plate_{p_idx+1}"
                for r_i, r_idx in enumerate(b):
                    row_let = row_letters[r_i]
                    for col_idx in range(5, df_survey_raw.shape[1]):
                        v_val = df_survey_raw.iloc[r_idx, col_idx]
                        if pd.notna(v_val):
                            survey_recs.append({
                                'Plate_1536': plt_name,
                                'Well_1536': f"{row_let}{(col_idx-4):02d}",
                                'Measured_Volume_uL': float(v_val)
                            })

            survey_df = pd.DataFrame(survey_recs)
            merged_pre = pd.merge(map_df, survey_df, on=['Plate_1536', 'Well_1536'], how='left')

            clean_wells = merged_pre[merged_pre['Measured_Volume_uL'] >= vol_cutoff].copy()
            depleted_wells = merged_pre[merged_pre['Measured_Volume_uL'] < vol_cutoff].copy()

            m1, m2, m3 = st.columns(3)
            m1.metric("Total Library Compounds", len(map_df))
            m2.metric("Sufficient Volume Wells (Passed)", len(clean_wells))
            m3.metric("Depleted Wells (Filtered Out)", len(depleted_wells))

            st.subheader("Filtered Master Map (Ready for Pooling Engine)")
            st.dataframe(clean_wells.head(15), use_container_width=True)

            buf_clean = io.StringIO()
            clean_wells.to_csv(buf_clean, index=False)
            st.download_button("⬇️ Download Cleaned 1536 Map (.csv)", buf_clean.getvalue(), "1536_master_map_sufficient_vol.csv", "text/csv", type="primary")

            if not depleted_wells.empty:
                st.subheader("⚠️ Depleted Compounds Reorder Manifest")
                st.dataframe(depleted_wells[['NCGC_ID', 'Plate_1536', 'Well_1536', 'Measured_Volume_uL']], use_container_width=True)
                
                buf_dep = io.StringIO()
                depleted_wells.to_csv(buf_dep, index=False)
                st.download_button("⬇️ Download Reorder Manifest (.csv)", buf_dep.getvalue(), "depleted_compounds_reorder_list.csv", "text/csv")

        except Exception as ex:
            st.error(f"Error parsing survey file: {ex}")

# ==========================================
# TAB 4: POST-RUN ECHO EXCEPTION RECONCILER
# ==========================================
with tab4:
    st.subheader("⚡ Post-Run Echo Exception Reconciler")
    st.write("Upload an Echo Exception/Transfer Report after a run to automatically calculate corrected DMSO back-flushes and strip skipped compounds from downstream manifests.")

    col_r1, col_r2 = st.columns(2)
    with col_r1:
        manifest_up = st.file_uploader("1. Upload Source Prep Manifest (.xlsx or .csv)", type=["xlsx", "csv"], key="m_up")
    with col_r2:
        report_up = st.file_uploader("2. Upload Echo Exception Report (.csv or .xlsx)", type=["csv", "xlsx"], key="r_up")

    target_vol_ul_recon = st.number_input("Target 384 Well Working Volume (µL)", min_value=2.0, max_value=50.0, value=10.0, step=1.0)
    aliquot_vol_nl_recon = st.number_input("Aliquot Volume per Compound (nL)", min_value=10, max_value=5000, value=500, step=100)

    if manifest_up is not None and report_up is not None:
        try:
            if manifest_up.name.endswith('.csv'): orig_df = pd.read_csv(manifest_up)
            else: orig_df = pd.read_excel(manifest_up)

            if report_up.name.endswith('.csv'): exc_df = pd.read_csv(report_up)
            else: exc_df = pd.read_excel(report_up)

            exc_df.columns = [str(c).strip() for c in exc_df.columns]
            dest_well_col = next((c for c in exc_df.columns if 'DEST' in c.upper() and 'WELL' in c.upper()), None)
            
            if dest_well_col:
                failed_counts = exc_df[exc_df[dest_well_col].notnull()].groupby(dest_well_col).size().to_dict()
                
                reconciled_rows = []
                for well_384, group in orig_df.groupby('Source_Well_384'):
                    assigned_cnt = len(group)
                    failed_cnt = failed_counts.get(well_384, 0)
                    actual_cnt = max(0, assigned_cnt - failed_cnt)
                    
                    actual_comp_fluid_nl = actual_cnt * aliquot_vol_nl_recon
                    target_total_nl = target_vol_ul_recon * 1000.0
                    corrected_backflush_nl = max(0.0, target_total_nl - actual_comp_fluid_nl)
                    
                    reconciled_rows.append({
                        'Source_Plate_384': group['Source_Plate_384'].iloc[0],
                        'Source_Well_384': well_384,
                        'Assigned_Compounds': assigned_cnt,
                        'Failed_Transfers': failed_cnt,
                        'Actual_Compounds_Added': actual_cnt,
                        'Actual_Compound_Fluid_nL': actual_comp_fluid_nl,
                        'Target_Total_Vol_nL': target_total_nl,
                        'Corrected_DMSO_Backflush_nL': corrected_backflush_nl
                    })

                recon_summary_df = pd.DataFrame(reconciled_rows)
                skewed_wells = recon_summary_df[recon_summary_df['Failed_Transfers'] > 0]

                rm1, rm2, rm3 = st.columns(3)
                rm1.metric("Total 384 Pools Evaluated", len(recon_summary_df))
                rm2.metric("Total Failed Echo Transfers", len(exc_df))
                rm3.metric("384 Wells Requiring Corrected Back-flush", len(skewed_wells))

                st.subheader("Corrected DMSO Back-flush Picklist (For Echo Top-Up Run)")
                backflush_picklist = skewed_wells[['Source_Plate_384', 'Source_Well_384', 'Corrected_DMSO_Backflush_nL']].copy()
                backflush_picklist.columns = ['Destination Plate Name', 'Destination Well', 'Transfer Volume']
                backflush_picklist['Source Plate Name'] = 'DMSO_SOURCE'
                backflush_picklist['Source Well'] = 'A01'
                backflush_picklist = backflush_picklist[['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume']]

                st.dataframe(backflush_picklist, use_container_width=True)

                buf_bf = io.StringIO()
                backflush_picklist.to_csv(buf_bf, index=False)
                st.download_button("⬇️ Download Corrected Back-flush Picklist (.csv)", buf_bf.getvalue(), "corrected_dmso_backflush_picklist.csv", "text/csv", type="primary")

        except Exception as ex_recon:
            st.error(f"Error reconciling report: {ex_recon}")
