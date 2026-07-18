import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import os
import io
import math
import json

st.set_page_config(page_title="ASMS Pooling Engine", page_icon="⚜️", layout="wide")

st.title("NCATS ASMS Compound Pooling Engine")
st.markdown("Created by Colin Gordy for use in the development of a semi-automated, small-molecule binders discovery assay utilizing HRMS. This tool automates multi-stage HTS workflows: Compiles 1536-well library entries into consolidated 384-well acoustic source pools, tracks volume normalization, and maps subsequent nanoliter transfers to 96-well assay target plates. Using .SDF files containing NCGC IDs and SMILES, the tool returns three outputs: (1) A .xlsx Echo script for creating an acoustic 384-well source plate with the 1536-well CoMa library plates, (2) A .xlsx Echo script for using the 384-well acoustic source plate and standard, 96-well KingFisher Flex plates as the destination, (3) A .HTML file for an interactive visualization of both plate maps. Parameters can be customized using the side toolbar.")
# ==========================================
# 1. Sidebar Control Panel
# ==========================================
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

# Main File Upload Area
st.markdown("### Upload Core Campaign Assets")
up_col1, up_col2 = st.columns(2)

with up_col1:
    uploaded_file = st.file_uploader("Required: Choose an SDF Library File", type=["sdf"])
    plate_prefix = st.text_input("Plate Name Prefix", value="ASMS")

with up_col2:
    # 🧪 NEW: Secondary, completely optional uploader for compound management inventories
    uploaded_inventory = st.file_uploader("Optional Upstream Link: Upload 1536 Master Plate Maps", type=["csv", "xlsx"], help="Provide the manifest file containing real-world freezer locations to generate the initial 1536 to 384 pool picklist file.")

# ==========================================
# 2. Computational Core Functions
# ==========================================
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
        <h1>NCATS ASMS Interactive Navigator</h1>
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
                let info = document.createElement('div'); info.className = 'compound-info';
                info.innerHTML = `
                    <div class="compound-id">${c.id}</div>
                    <div><strong>Exact Mass:</strong> ${c.mass.toFixed(4)} Da</div>
                    <div><strong>Target M/Z:</strong> ${c.mz.toFixed(4)}</div>
                    <div style="margin-top:5px; color:#64748b; font-size:11px; word-break:break-all;"><strong>SMILES:</strong> ${c.smiles}</div>
                `;
                let imgDiv = document.createElement('div'); imgDiv.className = 'struct-img';
                if(c.img) imgDiv.innerHTML = c.img;
                card.appendChild(info); card.appendChild(imgDiv); listContainer.appendChild(card);
            });
        }

        setViewType('source');
    </script>
</body>
</html>
"""
    return html_template.replace("{js_data_payload}", full_payload)

# ==========================================
# 3. Main System Pipeline Processing
# ==========================================
if uploaded_file is not None:
    max_possible_compound_vol_nl = pool_size * vol_per_comp
    target_source_vol_nl = target_source_vol_ul * 1000.0
    
    if max_possible_compound_vol_nl > target_source_vol_nl:
        st.error(f"❌ **Physical Fluidic Paradox Error:** You have requested a pool size of **{pool_size} compounds** at **{vol_per_comp} nL** each. This requires a minimum fluid volume of **{max_possible_compound_vol_nl / 1000.0} µL** per well from the library compound aliquots alone, which physically overflows your target Echo source well capacity of **{target_source_vol_ul} µL**. The liquid handler cannot execute negative DMSO back-flushes. Please reduce the compounds per well, lower the individual aliquot volume, or expand the target source well volume configuration.")
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
                well_idx = idx % 96
                assigned_96_well = assay_coordinates[well_idx]
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
            dash_col2.metric("384-Well Source Pools (Out 1)", total_384_wells)
            dash_col3.metric("96-Well Plates Needed", plates_needed)
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
            # 4. Multi-Deliverable Export Hub
            # ==========================================
            st.markdown("### Download Campaign Assets")
            
            # Dynamically adjusts layout grid columns depending on whether inventory maps are active
            if uploaded_inventory is not None:
                down_col0, down_col1, down_col2, down_col3, down_col4 = st.columns(5)
            else:
                down_col1, down_col2, down_col3, down_col4 = st.columns(4)
            
            # 🧪 NEW: Upstream 1536 -> 384 Echo consolidation processing pipeline
            if uploaded_inventory is not None:
                try:
                    if uploaded_inventory.name.endswith('.csv'):
                        inv_df = pd.read_csv(uploaded_inventory)
                    else:
                        inv_df = pd.read_excel(uploaded_inventory)
                    
                    # Cleans string headers for cross-reference matching
                    inv_df.columns = [str(c).strip() for c in inv_df.columns]
                    
                    # ⚠️ Placeholder matching keys: Assumes 'NCGC_ID', 'Plate_1536', 'Well_1536' column topology
                    expected_id_col = next((c for c in inv_df.columns if c.upper() in ['NCGC_ID', 'SAMPLE_ID', 'ID', 'COMPOUND_ID']), None)
                    expected_plt_col = next((c for c in inv_df.columns if 'PLATE' in c.upper() or 'BARCODE' in c.upper()), None)
                    expected_wel_col = next((c for c in inv_df.columns if 'WELL' in c.upper() or 'COORD' in c.upper()), None)
                    
                    if expected_id_col and expected_plt_col and expected_wel_col:
                        inv_df['match_id'] = inv_df[expected_id_col].astype(str).str.strip()
                        source_map['match_id'] = source_map['NCGC_ID'].astype(str).str.strip()
                        
                        consolidation_df = pd.merge(source_map, inv_df, on='match_id', how='inner')
                        
                        if not consolidation_df.empty:
                            picklist_1536_to_384 = consolidation_df[[
                                expected_plt_col, expected_wel_col, 'Source_Plate_384', 'Source_Well_384'
                            ]].copy()
                            picklist_1536_to_384['Transfer Volume'] = vol_per_comp
                            picklist_1536_to_384.columns = [
                                'Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well', 'Transfer Volume'
                            ]
                            
                            # Removes duplicates to create flat transfer matrix commands
                            picklist_1536_to_384 = picklist_1536_to_384.drop_duplicates().reset_index(drop=True)
                            
                            buf_up = io.StringIO()
                            picklist_1536_to_384.to_csv(buf_up, index=False)
                            csv_up = buf_up.getvalue()
                            
                            down_col0.download_button(
                                label="0. Download 1536 ➔ 384 Picklist (.csv)",
                                data=csv_up,
                                file_name=f"{plate_prefix.strip().lower()}_1536_to_384_consolidation.csv",
                                mime="text/csv",
                                use_container_width=True
                            )
                        else:
                            down_col0.error("0 Match IDs identified between files.")
                    else:
                        down_col0.error("Missing mapping column keys.")
                except Exception as ex:
                    down_col0.error(f"Upstream pipeline error: {str(ex)}")

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
