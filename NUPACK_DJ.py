import streamlit as st
import RNA
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
import streamlit.components.v1 as components
import os

# 🌟 尝试导入 NUPACK
try:
    from nupack import *
    nupack_available = True
except ImportError:
    nupack_available = False

st.set_page_config(page_title="核酸多链分析平台", layout="wide")

st.title("🧬 高级核酸序列分析与多链配对平台")
st.markdown("集成 ViennaRNA 与 NUPACK 算法，支持单链分析与多链杂交模拟。")

tab1, tab2 = st.tabs(["单链常规分析 (ViennaRNA)", "多链杂交模拟 (NUPACK)"])

# ==========================================
# 标签页 1：单链分析
# ==========================================
with tab1:
    st.subheader("模式一：单链二级结构与参数计算")
    sequence_input = st.text_area("输入单条核酸序列:", value="GCGCUUCGCCGCGCCCGUGCUG", height=100)
    
    if st.button("开始单链分析"):
        clean_seq = sequence_input.upper().replace(" ", "").replace("\n", "")
        if clean_seq:
            # 计算数据
            g_count = clean_seq.count('G')
            c_count = clean_seq.count('C')
            gc_content = (g_count + c_count) / len(clean_seq) * 100
            struct, mfe = RNA.fold(clean_seq)
            
            try:
                tm_value = mt.Tm_NN(clean_seq, nn_table=mt.RNA_NN1)
                tm_display = f"{tm_value:.2f} °C"
            except:
                tm_display = "无法计算"

            # 展示数据
            st.success("分析完成！")
            col_a, col_b, col_c = st.columns(3)
            col_a.metric("GC 含量", f"{gc_content:.2f}%")
            col_b.metric("最小自由能 (MFE)", f"{mfe:.2f} kcal/mol")
            col_c.metric("预测 Tm 值", tm_display)
            
            st.text_input("二级结构 (点号-括号):", value=struct)
            
            # 绘图
            temp_svg = "temp_single.svg"
            RNA.svg_rna_plot(clean_seq, struct, temp_svg)
            with open(temp_svg, "r") as f:
                svg_code = f.read()
            components.html(f"<div style='text-align: center;'>{svg_code}</div>", height=500)
            os.remove(temp_svg)

# ==========================================
# 标签页 2：NUPACK 多链分析
# ==========================================
with tab2:
    st.subheader("模式二：多链杂交模拟")
    if not nupack_available:
        st.warning("⚠️ 检测到当前环境未安装 NUPACK。部署到云端后将自动激活。")
    
    col_p1, col_p2, col_p3 = st.columns(3)
    with col_p1:
        n_mat = st.selectbox("材质:", ["RNA", "DNA"])
        n_temp = st.number_input("温度 (°C):", value=37.0)
    with col_p2:
        n_na = st.number_input("Na+ (M):", value=1.0)
        n_mg = st.number_input("Mg++ (M):", value=0.0)
    with col_p3:
        max_size = st.number_input("最大尺寸 (聚体):", min_value=1, max_value=4, value=2)

    st.markdown("####  反应组分与浓度")
    default_data = pd.DataFrame({
        "名称": ["链A", "链B"],
        "序列": ["AGUCUAGGAUUCGGCGUG", "CACGCCGAAUCCUAGACU"],
        "浓度 (µM)": [1.0, 1.0]
    })
    edited_df = st.data_editor(default_data, num_rows="dynamic", use_container_width=True)
    
    if st.button("启动 NUPACK 分析"):
        if not nupack_available:
            st.error("本地环境无法运行 NUPACK。")
        else:
            with st.spinner("计算中..."):
                try:
                    my_model = Model(material=n_mat, celsius=n_temp, sodium=n_na, magnesium=n_mg)
                    strands_dict = {}
                    
                    # 🌟 核心修复 1：创建一个我们自己的“备忘录字典”
                    seq_map = {} 
                    
                    for _, row in edited_df.iterrows():
                        s_name, s_seq = str(row["名称"]).strip(), str(row["序列"]).upper().replace(" ", "")
                        if s_seq and s_name:
                            strand_obj = Strand(s_seq, name=s_name)
                            strands_dict[strand_obj] = float(row["浓度 (µM)"]) * 1e-6
                            # 把名字和序列存进我们的备忘录
                            seq_map[s_name] = s_seq 

                    if not strands_dict:
                        st.warning("请输入有效序列！")
                        st.stop()

                    my_tube = Tube(strands=strands_dict, complexes=SetSpec(max_size=max_size), name="Tube1")
                    res = tube_analysis(tubes=[my_tube], model=my_model)[my_tube]
                    
                    results = []
                    # 遍历所有生成的复合物
                    for complex_item, conc in res.complex_concentrations.items():
                        conc_uM = conc / 1e-6
                        if conc_uM > 1e-4:
                            # 计算该复合物的 MFE 结构用于绘图
                            mfe_res = mfe(strands=complex_item.strands, model=my_model)[0]
                            results.append({
                                "复合物": complex_item.name,
                                "浓度 (µM)": conc_uM,
                                "MFE": mfe_res.energy,
                                "结构": str(mfe_res.structure),
                                "obj": complex_item, # 存储对象用于后续取名字
                                "struct_v": str(mfe_res.structure).replace("+", "&") # 转换符号
                            })
                    
                    df_res = pd.DataFrame(results).sort_values("浓度 (µM)", ascending=False).reset_index(drop=True)
                    st.success("计算完成！")
                    st.bar_chart(df_res.set_index("复合物")["浓度 (µM)"])
                    st.dataframe(df_res[["复合物", "浓度 (µM)", "MFE", "结构"]], use_container_width=True)

                    # 🌟 重点绘图区：按概率排列生成多个结构图
                    if not df_res.empty:
                        st.markdown("---")
                        st.subheader("🖼️ 产物结构分布图 (按生成概率排序)")
                        
                        # 设置想要展示的产物数量，建议前 3-5 个，避免页面过长
                        top_n = st.slider("展示排名前几位的产物图？", 1, 10, 3)
                        
                        # 遍历 DataFrame 的前 N 行
                        for i, row in df_res.head(top_n).iterrows():
                            # 计算百分比概率 (该复合物浓度 / 总浓度)
                            # 注意：这里简化为直接展示浓度，因为 NUPACK 的浓度直接反映了概率
                            prob_text = f"生成浓度: {row['平衡浓度 (µM)']:.4f} µM"
                            
                            with st.expander(f"排行 #{i+1}: {row['复合物']} ({prob_text})", expanded=(i==0)):
                                col_text, col_plot = st.columns([1, 2])
                                
                                with col_text:
                                    st.write(f"**能量 (MFE):** {row['MFE']:.2f} kcal/mol")
                                    st.write("**结构代码:**")
                                    st.code(row['结构'], language="text")
                                
                                with col_plot:
                                    # 拼接序列和结构（使用备忘录模式获取序列）
                                    combined_seq = "&".join([seq_map[s.name] for s in row["obj"].strands])
                                    vienna_struct = row["struct_v"]
                                    
                                    # 为每个图创建一个临时文件
                                    plot_file = f"temp_multi_{i}.svg"
                                    RNA.svg_rna_plot(combined_seq, vienna_struct, plot_file)
                                    
                                    with open(plot_file, "r") as f:
                                        svg_content = f.read()
                                    
                                    # 渲染 SVG
                                    st.components.v1.html(
                                        f"<div style='text-align:center;'>{svg_content}</div>", 
                                        height=450
                                    )
                                    
                                    # 清理文件
                                    if os.path.exists(plot_file):
                                        os.remove(plot_file)

                except Exception as e:
                    st.error(f"计算出错: {e}")