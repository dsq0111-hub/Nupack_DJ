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
st.markdown("集成 ViennaRNA 与 NUPACK (加州理工) 底层算法。")

# 使用标签页区分不同功能
tab1, tab2 = st.tabs(["单链常规分析 (ViennaRNA)", "多链复杂杂交模拟 (NUPACK)"])

# ==========================================
# 标签页 1：单链分析 (你之前的代码)
# ==========================================
with tab1:
    st.subheader("模式一：单链二级结构与参数计算")
    sequence_input = st.text_area("输入单条 RNA 序列:", value="GCGCUUCGCCGCGCCCGUGCUG", height=100)
    
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
with tab2:
    st.subheader("模式二：多条核酸链的试管级平衡态模拟 (NUPACK)")
    st.markdown("输入多条链,NUPACK 将模拟它们在试管中结合的最稳定状态。")
    
    if not nupack_available:
        st.warning("⚠️ 检测到当前环境未安装 NUPACK。由于 NUPACK 仅支持 Linux/Mac,此功能将在您部署到 Streamlit Cloud 后自动激活！")
    
    col1, col2 = st.columns(2)
    with col1:
        nupack_material = st.selectbox("核酸材质:", ["rna", "dna"])
        nupack_temp = st.number_input("模拟温度 (°C):", value=37.0)
    with col2:
        # 让用户输入多条链，每行一条
        multi_strands_input = st.text_area("输入多条序列（每行输入一条链）:", value="UGAGCUUCUGAGAAAUUCAA\nUUGAAUUUCUCAGAAGCUCA", height=120)
    
    if st.button("🔬 启动 NUPACK 多链模拟"):
        if not nupack_available:
            st.error("本地环境无法运行 NUPACK,请部署至云端后重试。")
        else:
            with st.spinner("NUPACK 引擎正在计算，请稍候..."):
                # 1. 获取并清理用户输入的每一条链
                raw_lines = multi_strands_input.strip().split("\n")
                clean_lines = [line.upper().replace(" ", "") for line in raw_lines if line.strip()]
                
                if len(clean_lines) < 2:
                    st.warning("请至少输入 2 条序列来进行多链分析！")
                else:
                    # 2. 设置 NUPACK 的热力学模型 (温度、材质)
                    my_model = Model(material=nupack_material, celsius=nupack_temp)
                    
                    # 3. 将输入的文本转化为 NUPACK 认识的 Strand 对象
                    strands_list = []
                    for i, seq in enumerate(clean_lines):
                        strands_list.append(Strand(seq, name=f"Strand_{i+1}"))
                    
                    # 4. 召唤 NUPACK 的 mfe 计算所有链结合的最稳定结构
                    # NUPACK 的 mfe 函数可以直接接受一个 Strand 列表
                    try:
                        mfe_results = mfe(strands=strands_list, model=my_model)
                        
                        # 获取能量最低（最可能）的那个结果
                        best_result = mfe_results[0]
                        
                        st.success("🎉 NUPACK 模拟完成！")
                        st.success("🎉 NUPACK 模拟完成！")
                        st.metric(label="复合体最稳定自由能 (MFE)", value=f"{best_result.energy:.2f} kcal/mol")
                        
                        # 1. 准备绘图所需的“翻译”数据
                        # NUPACK 的序列是一个列表，我们需要用 '&' 将它们连接起来
                        combined_seq = "&".join(clean_lines)
                        
                        # NUPACK 的结构字符串中使用 '+' 分隔链，我们需要换成 '&' 才能让 ViennaRNA 绘图
                        vienna_struct = str(best_result.structure).replace("+", "&")
                        
                        st.write("👉 最佳结合状态下的二维结构 (点号-括号表示):")
                        st.code(str(best_result.structure), language="text")
                        
                        # 2. 生成 SVG 矢量结构图
                        # 我们给多链图起个不同的名字，防止和单链图冲突
                        temp_multi_svg = "temp_multi_struct.svg"
                        
                        # 调用 ViennaRNA 的绘图函数
                        # 它非常聪明，看到序列和结构中有 '&' 就会自动按照多链杂交的模式画图
                        RNA.svg_rna_plot(combined_seq, vienna_struct, temp_multi_multi_svg)
                        
                        # 3. 在网页上渲染图片
                        st.subheader("🖼️ 多链杂交二级结构图")
                        with open(temp_multi_svg, "r") as f:
                            svg_code = f.read()
                        
                        # 居中展示
                        components.html(f"<div style='text-align: center;'>{svg_code}</div>", height=600)

                    except Exception as e:
                        st.error(f"NUPACK 计算出错：{e}")




            