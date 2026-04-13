import streamlit as st
import RNA
from Bio.SeqUtils import MeltingTemp as mt
import streamlit.components.v1 as components

# 🌟 尝试导入 NUPACK，如果是在本地 Windows 测试会报错，但在云端 Linux 会成功！
try:
    from nupack import *
    nupack_available = True
except ImportError:
    nupack_available = False

st.set_page_config(page_title="核酸分析平台", layout="wide")

st.title("🧬 高级核酸序列分析与多链配对平台")
st.markdown("集成 ViennaRNA 与 NUPACK 底层算法。")

# 使用标签页区分不同功能
tab1, tab2 = st.tabs(["单链常规分析 (ViennaRNA)", "多链复杂杂交模拟 (NUPACK)"])

# ==========================================
# 标签页 1：单链分析 (你之前的代码)
# ==========================================
with tab1:
    st.subheader("模式一：单链二级结构与参数计算")
    sequence_input = st.text_area("输入单条核酸序列:", value="GCGCUUCGCCGCGCCCGUGCUG", height=100)
    
    if st.button("开始单链分析"):
        clean_seq = sequence_input.upper().replace(" ", "").replace("\n", "")
        # ... 这里的代码就是上一节课的 GC含量、MFE、SVG 画图代码，你可以直接粘贴过来 ...   # ==========================================
        
        # (1) 计算 GC 含量
        g_count = clean_seq.count('G')
        c_count = clean_seq.count('C')
        gc_content = (g_count + c_count) / len(clean_seq) * 100
        
        # (2) 计算 MFE 和 点号-括号结构
        struct, mfe = RNA.fold(clean_seq)
        
        # (3) 计算 Tm 值
        # 注意：这里调用了 RNA 专属的最近邻热力学参数表 (mt.RNA_NN1)
        try:
            tm_value = mt.Tm_NN(clean_seq, nn_table=mt.RNA_NN1)
            tm_display = f"{tm_value:.2f} °C"
        except Exception as e:
            # 如果序列太短或有奇怪的字符，可能会报错，我们做个防护
            tm_display = "无法计算 (序列过短或含非标准碱基)"
            
        # (4) 生成 SVG 二级结构图，先临时保存在电脑里
        temp_svg_file = "temp_struct.svg"
        RNA.svg_rna_plot(clean_seq, struct, temp_svg_file)
        
        # ==========================================
        st.success("分析完成！")

        # 用三列并排显示三个核心数据
        col_a, col_b, col_c = st.columns(3)
        col_a.metric(label="GC 含量", value=f"{gc_content:.2f}%")
        col_b.metric(label="最小自由能 (MFE)", value=f"{mfe:.2f} kcal/mol")
        col_c.metric(label="预测 Tm 值", value=tm_display)
        
        # 显示序列结构（用 text_input 方便别人复制）
        st.text_input("二级结构 (点号-括号表示法):", value=struct)
        
        # 读取并在网页上渲染那张 SVG 图片
        st.subheader(" 预测二级结构图")
        with open(temp_svg_file, "r") as f:
            svg_code = f.read()
        components.html(f"<div style='text-align: center;'>{svg_code}</div>", height=500)
       

# ==========================================
# 标签页 2：NUPACK 多链杂交分析 (全新高级功能)
# ==========================================
# ==========================================
# 标签页 2：NUPACK 试管平衡态分析 (高度还原 NUPACK 官方功能)
# ==========================================
with tab2:
    st.subheader("模式二：多链试管平衡态分析 (NUPACK Tube Analysis)")
    st.markdown("还原 NUPACK 官方网站的试管模拟功能。支持任意数量的核酸链、自定义浓度、以及盐离子浓度调节。")
    
    if not nupack_available:
        st.warning("检测到当前环境未安装 NUPACK。部署到 Streamlit Cloud 后将自动激活。")
    
    # --- 1. NUPACK 核心物理参数设置区 ---
    st.markdown("#### 物理参数设置 (Physical Parameters)")
    col1, col2, col3 = st.columns(3)
    with col1:
        n_material = st.selectbox("核酸类型 (Material):", ["RNA", "DNA"])
        n_temp = st.number_input("温度 (Temperature °C):", value=37.0)
    with col2:
        # 盐离子浓度对核酸杂交极其重要，NUPACK 官网默认 Na+ 为 1.0M
        n_na = st.number_input("钠离子 Na+ 浓度 (M):", value=1.0, step=0.1)
        n_mg = st.number_input("镁离子 Mg++ 浓度 (M):", value=0.0, step=0.001, format="%.3f")
    with col3:
        # 允许的最大聚合物大小
        max_size = st.number_input("最大复合物尺寸 (Max Complex Size):", min_value=1, max_value=5, value=2, help="例如设为2，则最多计算两条链结合的双链；设为3则计算到三聚体。数值越大计算越慢。")

    # --- 2. 序列录入区 ---
    st.markdown("#### 序列与初始浓度 (Strands & Initial Concentrations)")
    st.info(" 你可以直接点击单元格修改序列，或者在表格底部添加/删除行。")
    
    import pandas as pd
    
    # 提供一个最基础的通用双链默认值
    default_strands = pd.DataFrame({
        "链名称 (Name)": ["Strand_A", "Strand_B"],
        "序列 (Sequence)": ["AGUCUAGGAUUCGGCGUG", "CACGCCGAAUCCUAGACU"],
        "初始浓度 (µM)": [1.0, 1.0] 
    })
    
    edited_df = st.data_editor(default_strands, num_rows="dynamic", use_container_width=True)
    
    # --- 3. 分析运行区 ---
    if st.button("启动 NUPACK 平衡态分析"):
        if not nupack_available:
            st.error("本地环境无法运行 NUPACK,请部署至云端后重试。")
        else:
            with st.spinner("NUPACK 引擎正在计算配分函数与平衡态浓度，请稍候..."):
                try:
                    # 将界面上收集到的物理参数传给 NUPACK Model
                    my_model = Model(
                        material=n_material, 
                        celsius=n_temp, 
                        sodium=n_na,       # 引入钠离子
                        magnesium=n_mg     # 引入镁离子
                    )
                    strands_dict = {}
                    
                    # 遍历表格获取用户输入
                    for index, row in edited_df.iterrows():
                        name = str(row["链名称 (Name)"]).strip()
                        seq = str(row["序列 (Sequence)"]).upper().replace(" ", "")
                        
                        if seq and name:
                            try:
                                conc_uM = float(row["初始浓度 (µM)"])
                                strand_obj = Strand(seq, name=name)
                                # NUPACK 计算使用 M (摩尔/升)，这里将 µM 换算为 M
                                strands_dict[strand_obj] = conc_uM * 1e-6 
                            except ValueError:
                                st.error(f"'{name}' 的浓度输入有误，必须为数字！")
                                st.stop()
                    
                    if len(strands_dict) == 0:
                        st.warning("请至少输入一条有效的序列！")
                    else:
                        # 创建试管进行分析
                        my_tube = Tube(strands=strands_dict, complexes=SetSpec(max_size=max_size), name="Simulation_Tube")
                        tube_results = tube_analysis(tubes=[my_tube], model=my_model)
                        res = tube_results[my_tube]
                        
                        results_data = []
                        for complex_item, complex_result in res.complexes.items():
                            conc_uM = complex_result.concentration / 1e-6 
                            
                            # 过滤掉极低浓度的微量副产物（比如低于 0.0001 µM）使结果更干净
                            if conc_uM > 1e-4:
                                comp_mfe = mfe(strands=complex_item.strands, model=my_model)[0]
                                results_data.append({
                                    "复合物组合": complex_item.name, 
                                    "平衡浓度 (µM)": conc_uM,
                                    "MFE 能量 (kcal/mol)": comp_mfe.energy,
                                    "最稳定二级结构": str(comp_mfe.structure)
                                })
                        
                        df_results = pd.DataFrame(results_data)
                        df_results = df_results.sort_values(by="平衡浓度 (µM)", ascending=False).reset_index(drop=True)
                        
                        st.success(" 分析完成！")
                        
                        # --- 4. 结果可视化输出 ---
                        st.markdown("平衡态最终浓度分布 (Equilibrium Concentrations)")
                        st.markdown("该图展示了在充分反应后，试管中各种游离单链、杂交双链或多聚体的最终浓度。")
                        chart_data = df_results.set_index("复合物组合")["平衡浓度 (µM)"]
                        st.bar_chart(chart_data)
                        
                        st.markdown("### 详细复合物热力学属性 (Complex Thermodynamics)")
                        st.dataframe(df_results, use_container_width=True)

                except Exception as e:
                    st.error(f"NUPACK 计算异常，请检查序列中是否包含非标碱基：{e}")



            