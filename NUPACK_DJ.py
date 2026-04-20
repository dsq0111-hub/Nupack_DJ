import streamlit as st
import RNA
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
import streamlit.components.v1 as components
import os
import re

# 🌟 尝试导入 NUPACK
try:
    from nupack import *
    nupack_available = True
except ImportError:
    nupack_available = False

st.set_page_config(page_title="高级核酸分析平台", layout="wide", page_icon="🧬")
st.title("🧬 高级核酸序列分析与多链配对平台")

tab1, tab2 = st.tabs(["单链常规分析 (ViennaRNA)", "多链杂交模拟 (NUPACK)"])

# ==========================================
# 标签页 1：单链分析 (回调自动排版)
# ==========================================
def format_single_seq():
    raw = st.session_state.seq1_raw
    clean = raw.upper().replace(" ", "").replace("\n", "")
    st.session_state.seq1_raw = " ".join([clean[i:i+6] for i in range(0, len(clean), 6)])

with tab1:
    st.subheader("模式一：单链二级结构与参数计算")
    
    if 'seq1_raw' not in st.session_state:
        st.session_state.seq1_raw = "GCGCUU CGCCGC CCCGUG CUG"
        
    st.info("💡 提示：在此处随意输入序列 (支持小写、带空格)。输入完成后点击框外，系统会自动排版为标准格式")
    
    sequence_input = st.text_area(
        "输入单条 RNA/DNA 序列:", 
        key="seq1_raw", 
        on_change=format_single_seq, 
        height=100
    )
    
    st.markdown("####  底层配对算法控制")
    strict_wc = st.checkbox(" 严格禁止 G-U / G-T 摆动配对 (强制仅允许 A-T/U, G-C 经典配对)", value=True)
    
    if st.button("开始单链分析"):
        clean_seq = sequence_input.replace(" ", "")
        if clean_seq:
            # 算法控制阀门
            if strict_wc:
                RNA.cvar.noGU = 1 
            else:
                RNA.cvar.noGU = 0 
                
            g_count = clean_seq.count('G')
            c_count = clean_seq.count('C')
            gc_content = (g_count + c_count) / len(clean_seq) * 100
            
            struct, mfe_val = RNA.fold(clean_seq)
            
            try:
                tm_value = mt.Tm_NN(clean_seq, nn_table=mt.RNA_NN1)
                tm_display = f"{tm_value:.2f} °C"
            except Exception:
                tm_display = "无法计算"

            st.success("分析完成！")
            col_a, col_b, col_c = st.columns(3)
            col_a.metric(" GC 含量", f"{gc_content:.2f}%")
            col_b.metric(" 最小自由能 (MFE)", f"{mfe_val:.2f} kcal/mol")
            col_c.metric(" 预测 Tm 值", tm_display)
            
            st.text_input(" 二级结构 (点号-括号):", value=struct)
            
            temp_svg = "temp_single.svg"
            RNA.svg_rna_plot(clean_seq, struct, temp_svg)
            with open(temp_svg, "r") as f:
                svg_code = f.read()
            st.components.v1.html(
                f"<div style='text-align: center; background-color: white; border-radius: 10px; padding: 10px;'>{svg_code}</div>", 
                height=500
            )
            st.download_button(
                label=" 下载单链结构矢量图 (.svg)",
                data=svg_code,
                file_name="Single_Strand_Structure.svg",
                mime="image/svg+xml",
                key="dl_btn_single"
            )
            if os.path.exists(temp_svg):
                os.remove(temp_svg)

# ==========================================
# 标签页 2：NUPACK 多链模拟 (彻底修复内存丢失 Bug 版)
# ==========================================
with tab2:
    st.subheader("模式二：多链杂交平衡态模拟")
    if not nupack_available:
        st.warning(" 检测到当前环境未安装 NUPACK。部署到云端后将自动激活。")

def polish_svg(svg_str, chain_sequences):
        # 1. 抹除画蛇添足的 '&' 符号
        svg_str = re.sub(r'<text[^>]*>&amp;</text>', '', svg_str)
        svg_str = re.sub(r'<text[^>]*>&</text>', '', svg_str)
        
        # 2. 为不同的链分配科研配色
        colors = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#7E6148"]
        chain_lengths = [len(seq) for seq in chain_sequences]
        tracker = {"chain_idx": 0, "base_count": 0}
        
        # 新增：用来存储 5' 和 3' 的 SVG 标签，等下统一注入
        end_labels = []
        
        def color_and_label_injector(match):
            tag = match.group(0) # 获取整个 <text ...>A</text> 标签
            curr_chain = tracker["chain_idx"]
            color = colors[curr_chain % len(colors)]
            
            # 提取当前碱基的 X 和 Y 坐标
            x_match = re.search(r'x="([-\d\.]+)"', tag)
            y_match = re.search(r'y="([-\d\.]+)"', tag)
            
            if x_match and y_match:
                x = float(x_match.group(1))
                y = float(y_match.group(1))
                
                # 如果是这条链的【第一个碱基】(5' 端)
                if tracker["base_count"] == 0:
                    # 在它左上角偏移一点的位置，画一个 5'
                    end_labels.append(f"<text x='{x-12}' y='{y-12}' fill='{color}' font-size='12' font-weight='bold' font-family='Arial'>5'</text>")
                
                # 如果是这条链的【最后一个碱基】(3' 端)
                if tracker["base_count"] == chain_lengths[curr_chain] - 1:
                    # 在它右下角偏移一点的位置，画一个 3'
                    end_labels.append(f"<text x='{x+12}' y='{y+12}' fill='{color}' font-size='12' font-weight='bold' font-family='Arial'>3'</text>")
            
            tracker["base_count"] += 1
            # 链切换逻辑
            if tracker["base_count"] >= chain_lengths[curr_chain]:
                tracker["chain_idx"] += 1
                tracker["base_count"] = 0
                
            return tag.replace("<text ", f"<text fill='{color}' font-weight='bold' ")
            
        # 第一步：替换碱基颜色，并收集 5'/3' 坐标
        modified_svg = re.sub(r'<text[^>]*>[A-Za-z]</text>', color_and_label_injector, svg_str)
        
        # 第二步：将收集到的 5' 和 3' 标记，在 SVG 结束前强行插入进去
        if end_labels:
            labels_str = "\n".join(end_labels)
            modified_svg = modified_svg.replace("</svg>", f"{labels_str}\n</svg>")
            
        return modified_svg


    # 🌟 核心修复 1：定义绝对唯一的“主数据源 (Master Data)”
    if 'master_df' not in st.session_state:
        st.session_state.master_df = pd.DataFrame({
            "名称": ["Target", "Probe"],
            "序列": ["AGUCUAGGAUUCGGCGUG", "CACGCCGAAUCCUAGACU"],
            "浓度 (µM)": [1.0, 1.0]
        })
    if 'editor_key' not in st.session_state:
        st.session_state.editor_key = 0
    if 'nupack_results' not in st.session_state:
        st.session_state.nupack_results = None
    if 'nupack_seq_map' not in st.session_state:
        st.session_state.nupack_seq_map = {}

    col_p1, col_p2, col_p3 = st.columns(3)
    with col_p1:
        n_mat = st.selectbox("材质:", ["RNA", "DNA"])
    with col_p2:
        n_temp = st.number_input("温度 (°C):", value=37.0)
    with col_p3:
        max_size = st.number_input("最大尺寸 (几聚体):", min_value=1, max_value=4, value=2)
        n_na = 1.0 
        n_mg = 0.0

    st.markdown("---")
    st.markdown("####  实验记录与序列管理")
    
    col_file1, col_file2 = st.columns(2)
    with col_file1:
        uploaded_file = st.file_uploader("📂 导入历史序列数据 (.csv)", type=["csv"])

    # 🌟 核心修复 2：稳健的文件上传覆盖逻辑
    if uploaded_file is not None:
        if 'last_uploaded' not in st.session_state or st.session_state.last_uploaded != uploaded_file.name:
            try:
                st.session_state.master_df = pd.read_csv(uploaded_file)
                st.session_state.last_uploaded = uploaded_file.name
                st.session_state.editor_key += 1 # 强制刷新表格以显示新文件内容
            except Exception:
                st.error("读取文件失败，请检查格式。")

    # 🌟 核心修复 3：绝对安全的格式化逻辑 (防越界、防空行)
    if st.button(" 一键排版表格序列 (去空 / 大写 / 6位分隔)"):
        df = st.session_state.master_df.copy()
        for idx in df.index: # 使用 index 而不是 range(len)，防止删行后报错
            val = df.at[idx, "序列"]
            if pd.isna(val): # 防止用户留空行引发 NaN 报错
                df.at[idx, "序列"] = ""
            else:
                seq = str(val).upper().replace(" ", "").replace("\n", "")
                df.at[idx, "序列"] = " ".join([seq[j:j+6] for j in range(0, len(seq), 6)])
        st.session_state.master_df = df
        st.session_state.editor_key += 1 # 强制刷新表格 UI

    st.markdown("####  反应组分与浓度")
    
    # 渲染带有动态 Key 的表格，并实时将结果写回 master_df
    edited_df = st.data_editor(
        st.session_state.master_df, 
        key=f"editor_{st.session_state.editor_key}", 
        num_rows="dynamic", 
        use_container_width=True
    )
    st.session_state.master_df = edited_df

    with col_file2:
        csv_data = edited_df.to_csv(index=False).encode('utf-8-sig')
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.download_button(" 保存当前表格为存档 (.csv)", data=csv_data, file_name="NUPACK_Sequences.csv", mime="text/csv")

    if st.button(" 启动 NUPACK 分析"):
        if not nupack_available:
            st.error("本地环境无法运行 NUPACK。")
        else:
            with st.spinner("NUPACK 计算中..."):
                try:
                    my_model = Model(material=n_mat, celsius=n_temp, sodium=n_na, magnesium=n_mg)
                    strands_dict = {}
                    local_seq_map = {} 
                    
                    for _, row in edited_df.iterrows():
                        s_name = str(row["名称"]).strip()
                        val = row["序列"]
                        # 防御性代码：忽略空行和无效数据
                        s_seq = "" if pd.isna(val) else str(val).upper().replace(" ", "")
                        if s_seq and s_name and s_name != 'nan':
                            strands_dict[Strand(s_seq, name=s_name)] = float(row["浓度 (µM)"]) * 1e-6
                            local_seq_map[s_name] = s_seq 

                    if strands_dict:
                        my_tube = Tube(strands=strands_dict, complexes=SetSpec(max_size=max_size), name="Tube1")
                        res = tube_analysis(tubes=[my_tube], model=my_model)[my_tube]
                        
                        results = []
                        for complex_item, conc in res.complex_concentrations.items():
                            conc_uM = conc / 1e-6
                            if conc_uM > 1e-4:
                                mfe_res = mfe(strands=complex_item.strands, model=my_model)[0]
                                results.append({
                                    "复合物": complex_item.name,
                                    "浓度 (µM)": conc_uM,
                                    "MFE": mfe_res.energy,
                                    "结构": str(mfe_res.structure),
                                    "obj": complex_item, 
                                    "struct_v": str(mfe_res.structure).replace("+", "&") 
                                })
                        
                        st.session_state.nupack_results = pd.DataFrame(results).sort_values("浓度 (µM)", ascending=False).reset_index(drop=True)
                        st.session_state.nupack_seq_map = local_seq_map
                        st.success("🎉 计算完成！")
                    else:
                        st.warning("表格中没有有效的序列，请补充后再试！")
                except Exception as e:
                    st.error(f"计算出错: {e}")

    if st.session_state.nupack_results is not None:
        df_res = st.session_state.nupack_results
        seq_map = st.session_state.nupack_seq_map
        
        st.markdown("---")
        st.subheader(" 试管平衡态：产物分布统计")

        total_conc = df_res["浓度 (µM)"].sum()
        df_res["形成概率 (%)"] = (df_res["浓度 (µM)"] / total_conc) * 100 if total_conc > 0 else 0
        display_df = df_res[["复合物", "形成概率 (%)", "浓度 (µM)", "MFE", "结构"]]

        st.dataframe(
            display_df,
            column_config={
                "形成概率 (%)": st.column_config.ProgressColumn("形成概率 (%)", format="%.2f%%", min_value=0, max_value=100),
                "浓度 (µM)": st.column_config.NumberColumn(format="%.4f µM")
            },
            use_container_width=True, hide_index=True
        )

        st.markdown("###  产物结构分布图 (高清彩色版)")
        max_items = len(df_res)
        
        if max_items == 0:
            st.info(" 当前体系未生成显著的稳定复合物结构。")
        else:
            top_n = 1 if max_items == 1 else st.slider("展示排名前几位的产物图？", 1, max_items, min(3, max_items))
            
            for i, row in df_res.head(top_n).iterrows():
                prob_text = f"生成概率: {row['形成概率 (%)']:.2f}% | 浓度: {row['浓度 (µM)']:.4f} µM"
                
                with st.expander(f"排行 #{i+1}: {row['复合物']} ({prob_text})", expanded=(i==0)):
                    col_text, col_plot = st.columns([1, 2])
                    
                    with col_text:
                        st.write(f"**能量 (MFE):** {row['MFE']:.2f} kcal/mol")
                        st.write("**结构代码:**")
                        st.code(row['结构'], language="text")
                        
                        st.markdown("**链颜色说明:**")
                        colors_list = [" 红", " 蓝", " 绿", " 紫", " 橙"]
                        for idx, s in enumerate(row["obj"].strands):
                            st.caption(f"{colors_list[idx % len(colors_list)]} : {s.name}")
                    
                    with col_plot:
                        chain_sequences = [seq_map[s.name] for s in row["obj"].strands]
                        combined_seq = "&".join(chain_sequences)
                        plot_file = f"temp_multi_{i}.svg"
                        
                        RNA.svg_rna_plot(combined_seq, row["struct_v"], plot_file)
                        with open(plot_file, "r") as f:
                            raw_svg = f.read()
                            
                        polished_svg = polish_svg(raw_svg, chain_sequences)
                        st.components.v1.html(
                            f"<div style='text-align:center; background-color: white; border-radius: 10px; padding: 10px;'>{polished_svg}</div>", 
                            height=400
                        )
                        st.download_button(
                            label=" 下载此结构矢量图 (.svg)",
                            data=polished_svg, file_name=f"Rank{i+1}_{row['复合物']}.svg", mime="image/svg+xml", key=f"dl_btn_{i}"
                        )
                        if os.path.exists(plot_file):
                            os.remove(plot_file)