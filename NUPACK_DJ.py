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

st.set_page_config(page_title="核酸分析平台", layout="wide", page_icon="🧬")
st.title("🧬核酸序列分析与多链配对平台")

tab1, tab2 = st.tabs(["单链常规分析 (ViennaRNA)", "多链杂交模拟 (NUPACK)"])

# ==========================================
# 自动排版回调函数：去空格 -> 转大写 -> 6个一组
# ==========================================
def format_single_seq():
    raw = st.session_state.seq1_raw
    clean = raw.upper().replace(" ", "").replace("\n", "")
    st.session_state.seq1_raw = " ".join([clean[i:i+6] for i in range(0, len(clean), 6)])



# ==========================================
# 标签页 1：单链分析 
# ==========================================
with tab1:
    st.subheader("模式一：单链二级结构与参数计算")
    
    # 绑定回调函数到 Session State
    if 'seq1_raw' not in st.session_state:
        st.session_state.seq1_raw = "GCGCUU CGCCGC CCCGUG CUG"
        
    st.info("💡 提示：在此处随意输入序列 (支持小写、带空格)。输入完成后点击框外，系统会自动排版为标准格式")
    
    # 输入框，绑定了 on_change 自动排版事件
    sequence_input = st.text_area(
        "输入单条 RNA/DNA 序列:", 
        key="seq1_raw", 
        on_change=format_single_seq, 
        height=100
    )
    
    # 🌟 核心改进：底层算法逻辑控制器
    st.markdown("#### ⚙️ 底层配对算法控制")
    strict_wc = st.checkbox("🚫 严格禁止 G-U / G-T 摆动配对 (强制仅允许 A-T/U, G-C 经典配对)", value=True)
    
    if st.button("开始单链分析"):
        clean_seq = sequence_input.replace(" ", "")
        if clean_seq:
            # 🌟 核心：修改 ViennaRNA 底层全局变量
            if strict_wc:
                RNA.cvar.noGU = 1  # 彻底封杀 G-U 配对
            else:
                RNA.cvar.noGU = 0  # 恢复物理默认允许
                
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
            col_a.metric("GC 含量", f"{gc_content:.2f}%")
            col_b.metric("最小自由能 (MFE)", f"{mfe_val:.2f} kcal/mol")
            col_c.metric("预测 Tm 值", tm_display)
            
            st.text_input("二级结构 (点号-括号):", value=struct)
            
            temp_svg = "temp_single.svg"
            RNA.svg_rna_plot(clean_seq, struct, temp_svg)
            with open(temp_svg, "r") as f:
                svg_code = f.read()
            st.components.v1.html(
                f"<div style='text-align: center; background-color: white; border-radius: 10px; padding: 10px;'>{svg_code}</div>", 
                height=500
            )
            st.download_button(
                label="下载单链结构矢量图 (.svg)",
                data=svg_code,
                file_name="Single_Strand_Structure.svg",
                mime="image/svg+xml",
                key="dl_btn_single"
            )
            if os.path.exists(temp_svg):
                os.remove(temp_svg)

# ==========================================
# 标签页 2：NUPACK 多链模拟 
# ==========================================
with tab2:
    st.subheader("模式二：多链杂交平衡态模拟")
    if not nupack_available:
        st.warning("⚠️ 检测到当前环境未安装 NUPACK。部署到云端后将自动激活。")

    def polish_svg(svg_str, chain_sequences):
        svg_str = re.sub(r'<text[^>]*>&amp;</text>', '', svg_str)
        svg_str = re.sub(r'<text[^>]*>&</text>', '', svg_str)
        colors = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"]
        chain_lengths = [len(seq) for seq in chain_sequences]
        tracker = {"chain_idx": 0, "base_count": 0}
        
        def color_injector(match):
            tag = match.group(0)
            curr_chain = tracker["chain_idx"]
            color = colors[curr_chain % len(colors)]
            tracker["base_count"] += 1
            if tracker["base_count"] >= chain_lengths[curr_chain]:
                tracker["chain_idx"] += 1
                tracker["base_count"] = 0
            return tag.replace("<text ", f"<text fill='{color}' font-weight='bold' ")
            
        return re.sub(r'<text[^>]*>[A-Za-z]</text>', color_injector, svg_str)

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
        n_na = 1.0  # 隐藏不常用的盐离子，保持界面清爽
        n_mg = 0.0

    st.markdown("---")
    st.markdown("#### 实验记录与序列管理")
    
    col_file1, col_file2 = st.columns(2)
    with col_file1:
        uploaded_file = st.file_uploader("导入历史序列数据 (.csv)", type=["csv"])

    # 绑定表格数据到 Session State 方便格式化
    if 'nupack_input_df' not in st.session_state:
        if uploaded_file is not None:
            st.session_state.nupack_input_df = pd.read_csv(uploaded_file)
        else:
            st.session_state.nupack_input_df = pd.DataFrame({
                "名称": ["Target", "Probe"],
                "序列": ["AGUCUA GGAUUC GGCGUG", "CACGCC GAAUCC UAGACU"],
                "浓度 (µM)": [1.0, 1.0]
            })

    # 🌟 新增：表格 UI 强制刷新计数器
    if 'editor_key' not in st.session_state:
        st.session_state.editor_key = 0

    if st.button("✨ 一键排版表格序列 (去空 / 大写 / 6位分隔)"):
        # 使用 copy() 隔绝内存冲突
        df = st.session_state.nupack_input_df.copy()
        for i in range(len(df)):
            seq = str(df.at[i, "序列"]).upper().replace(" ", "").replace("\n", "")
            df.at[i, "序列"] = " ".join([seq[j:j+6] for j in range(0, len(seq), 6)])
        st.session_state.nupack_input_df = df
        # 核心：给计数器加 1，强制换一个新的表格组件
        st.session_state.editor_key += 1 

    st.markdown("#### 反应组分与浓度")
    
    # 🌟 修复 Bug：去掉了容易崩溃的旧 key，换成了带计数器的动态 key
    edited_df = st.data_editor(
        st.session_state.nupack_input_df, 
        key=f"editor_{st.session_state.editor_key}", 
        num_rows="dynamic", 
        use_container_width=True
    )
    
    # 将用户在网页上对表格做的任何手写修改，实时存入大脑
    st.session_state.nupack_input_df = edited_df


    with col_file2:
        csv_data = edited_df.to_csv(index=False).encode('utf-8-sig')
        st.markdown("<br><br>", unsafe_allow_html=True)
        st.download_button("保存当前表格为存档 (.csv)", data=csv_data, file_name="NUPACK_Sequences.csv", mime="text/csv")

    if st.button("🚀 启动 NUPACK 分析"):
        if not nupack_available:
            st.error("本地环境无法运行 NUPACK。")
        else:
            with st.spinner("NUPACK 计算中..."):
                try:
                    my_model = Model(material=n_mat, celsius=n_temp, sodium=n_na, magnesium=n_mg)
                    strands_dict = {}
                    local_seq_map = {} 
                    
                    for _, row in edited_df.iterrows():
                        s_name, s_seq = str(row["名称"]).strip(), str(row["序列"]).upper().replace(" ", "")
                        if s_seq and s_name:
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
                except Exception as e:
                    st.error(f"计算出错: {e}")

    if st.session_state.nupack_results is not None:
        df_res = st.session_state.nupack_results
        seq_map = st.session_state.nupack_seq_map
        
        st.markdown("---")
        st.subheader("产物分布统计")

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

        st.markdown("###产物结构分布图")
        max_items = len(df_res)
        
        if max_items == 0:
            st.info("⚠️ 当前体系未生成显著的稳定复合物结构。")
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
                        colors_list = ["🔴 红", "🔵 蓝", "🟢 绿", "🟣 紫", "🟠 橙"]
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
                            label="下载此结构矢量图 (.svg)",
                            data=polished_svg, file_name=f"Rank{i+1}_{row['复合物']}.svg", mime="image/svg+xml", key=f"dl_btn_{i}"
                        )
                        if os.path.exists(plot_file):
                            os.remove(plot_file)