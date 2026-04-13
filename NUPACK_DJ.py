import streamlit as st
import RNA
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
import streamlit.components.v1 as components
import os
import re  # 用于美化 SVG 图片

# 🌟 尝试导入 NUPACK
try:
    from nupack import *
    nupack_available = True
except ImportError:
    nupack_available = False

# ==========================================
# 页面全局设置
# ==========================================
st.set_page_config(page_title="高级核酸分析平台", layout="wide", page_icon="🧬")

st.title("🧬 高级核酸序列分析与多链配对平台")
st.markdown("集成 ViennaRNA 与 NUPACK核心算法。")

# 使用标签页区分不同功能
tab1, tab2 = st.tabs(["单链常规分析 (ViennaRNA)", "多链杂交模拟 (NUPACK)"])

# ==========================================
# 标签页 1：单链分析 (保持稳定)
# ==========================================
with tab1:
    st.subheader("模式一：单链二级结构与参数计算")
    sequence_input = st.text_area("输入单条 RNA/DNA 序列:", value="GCGCUUCGCCGCGCCCGUGCUG", height=100)
    
    if st.button("开始单链分析"):
        clean_seq = sequence_input.upper().replace(" ", "").replace("\n", "")
        if clean_seq:
            # 计算基础数据
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
        
            # 画图并展示 
            # ----------------------------------------
            temp_svg = "temp_single.svg"
            RNA.svg_rna_plot(clean_seq, struct, temp_svg)
            
            with open(temp_svg, "r") as f:
                svg_code = f.read()
            
            # 1. 在网页上居中渲染展示图片 (加了白色背景，看着更干净)
            st.components.v1.html(
                f"<div style='text-align: center; background-color: white; border-radius: 10px; padding: 10px;'>{svg_code}</div>", 
                height=500
            )
            
            # 2. 👇 新增的下载按钮 👇
            st.download_button(
                label=" 下载单链结构矢量图 (.svg)",
                data=svg_code,
                file_name="Single_Strand_Structure.svg",
                mime="image/svg+xml",
                key="dl_btn_single"
            )
            
            # 3. 清理临时文件
            if os.path.exists(temp_svg):
                os.remove(temp_svg)

# ==========================================
# 标签页 2：NUPACK 多链杂交
# ==========================================
with tab2:
    st.subheader("模式二：多链杂交模拟")
    if not nupack_available:
        st.warning("⚠️ 检测到当前环境未安装 NUPACK。部署到云端后将自动激活。")

    # ----------------------------------------
    # 辅助魔法函数：美化 SVG
    # ----------------------------------------
    def polish_svg(svg_str, chain_sequences):
        # 1. 抹除画蛇添足的 '&' 符号
        svg_str = re.sub(r'<text[^>]*>&amp;</text>', '', svg_str)
        svg_str = re.sub(r'<text[^>]*>&</text>', '', svg_str)
        
        # 2. 为不同的链分配科研配色
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

    # ----------------------------------------
    # 初始化缓存 (极度重要：解决滑块刷新问题)
    # ----------------------------------------
    if 'nupack_results' not in st.session_state:
        st.session_state.nupack_results = None
    if 'nupack_seq_map' not in st.session_state:
        st.session_state.nupack_seq_map = {}

    # ----------------------------------------
    # 参数与数据输入区
    # ----------------------------------------
    col_p1, col_p2, col_p3 = st.columns(3)
    with col_p1:
        n_mat = st.selectbox("材质:", ["RNA", "DNA"])
        n_temp = st.number_input("温度 (°C):", value=37.0)
    with col_p2:
        n_na = st.number_input("Na+ (M):", value=1.0)
        n_mg = st.number_input("Mg++ (M):", value=0.0)
    with col_p3:
        max_size = st.number_input("最大尺寸 (几聚体):", min_value=1, max_value=4, value=2)

    st.markdown("#### 反应组分与浓度")
    default_data = pd.DataFrame({
        "名称": ["链1", "链2"],
        "序列": ["AGUCUAGGAUUCGGCGUG", "CACGCCGAAUCCUAGACU"],
        "浓度 (µM)": [1.0, 1.0]
    })
    edited_df = st.data_editor(default_data, num_rows="dynamic", use_container_width=True)
    
    # ----------------------------------------
    # 逻辑核心：点击按钮只负责"计算"并"存入缓存"
    # ----------------------------------------
    if st.button("启动 NUPACK 分析"):
        if not nupack_available:
            st.error("本地环境无法运行 NUPACK。")
        else:
            with st.spinner("NUPACK 底层引擎计算中，请稍候..."):
                try:
                    my_model = Model(material=n_mat, celsius=n_temp, sodium=n_na, magnesium=n_mg)
                    strands_dict = {}
                    local_seq_map = {} 
                    
                    for _, row in edited_df.iterrows():
                        s_name, s_seq = str(row["名称"]).strip(), str(row["序列"]).upper().replace(" ", "")
                        if s_seq and s_name:
                            strands_dict[Strand(s_seq, name=s_name)] = float(row["浓度 (µM)"]) * 1e-6
                            local_seq_map[s_name] = s_seq 

                    if not strands_dict:
                        st.warning("请输入有效序列！")
                        st.stop()

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
                    
                    # 保存到缓存区
                    st.session_state.nupack_results = pd.DataFrame(results).sort_values("浓度 (µM)", ascending=False).reset_index(drop=True)
                    st.session_state.nupack_seq_map = local_seq_map
                    st.success("🎉 计算完成！数据已存入缓存。向下滚动查看结果。")

                except Exception as e:
                    st.error(f"计算出错: {e}")

    # ----------------------------------------
    # 渲染展示区：从缓存中读取数据进行画图
    # (这部分在按钮之外，所以拖动滑块时它不会消失)
    # ----------------------------------------
    if st.session_state.nupack_results is not None:
        df_res = st.session_state.nupack_results
        seq_map = st.session_state.nupack_seq_map
        
        st.markdown("---")
        st.subheader("试管平衡态数据")
        st.bar_chart(df_res.set_index("复合物")["浓度 (µM)"])
        st.dataframe(df_res[["复合物", "浓度 (µM)", "MFE", "结构"]], use_container_width=True)

        st.markdown("###产物结构分布图 ")
        
        # 滑块控件：限制最大范围为产物总数
        max_items = len(df_res)
        default_val = min(3, max_items)
        top_n = st.slider("展示排名前几位的产物图？", 1, max_items, default_val)
        
        for i, row in df_res.head(top_n).iterrows():
            prob_text = f"生成浓度: {row['浓度 (µM)']:.4f} µM"
            
            with st.expander(f"排行 #{i+1}: {row['复合物']} ({prob_text})", expanded=(i==0)):
                col_text, col_plot = st.columns([1, 2])
                
                with col_text:
                    st.write(f"**能量 (MFE):** {row['MFE']:.2f} kcal/mol")
                    st.write("**结构代码:**")
                    st.code(row['结构'], language="text")
                    
                    # 图例说明
                    st.markdown("**颜色说明:**")
                    colors_list = ["🔴 红", "🔵 蓝", "🟢 绿", "🟣 紫", "🟠 橙"]
                    for idx, s in enumerate(row["obj"].strands):
                        st.caption(f"{colors_list[idx % len(colors_list)]} : {s.name}")
                
                with col_plot:
                    chain_sequences = [seq_map[s.name] for s in row["obj"].strands]
                    combined_seq = "&".join(chain_sequences)
                    vienna_struct = row["struct_v"]
                    plot_file = f"temp_multi_{i}.svg"
                    
                    # 画原始图
                    RNA.svg_rna_plot(combined_seq, vienna_struct, plot_file)
                    with open(plot_file, "r") as f:
                        raw_svg = f.read()
                        
                    # 美化图片
                    polished_svg = polish_svg(raw_svg, chain_sequences)
                    
                    # 展示图片
                    st.components.v1.html(
                        f"<div style='text-align:center; background-color: white; border-radius: 10px; padding: 10px;'>{polished_svg}</div>", 
                        height=400
                    )
                    
                    # 提供下载按钮
                    st.download_button(
                        label="下载此结构矢量图 (.svg)",
                        data=polished_svg,
                        file_name=f"Rank{i+1}_{row['复合物']}_Structure.svg",
                        mime="image/svg+xml",
                        key=f"dl_btn_{i}"
                    )
                    
                    if os.path.exists(plot_file):
                        os.remove(plot_file)