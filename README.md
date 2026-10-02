# task26100305

含重复值的 Mann-Whitney 精确检验与位置差（A − B）置信集合反演，纯标准库实现
（Python ≥ 3.9，仅用 `fractions` / `itertools`），所有计算均为精确有理数运算。

## 接口

```python
from mannwhitney import exact_test, invert_shift

# 精确检验：输入为两组非空十进制字符串，按精确数值解释
r = exact_test(["1.5", "2", "2", "3.25"], ["2", "2.5"])
r.u        # Fraction(3, 1)   —— 观测 U（A>B 计 1，相等计 1/2）
r.p_value  # Fraction(14, 15) —— 精确双侧 p 值（全标签重枚举）

# 置信集合反演：置信度支持 0.90 / 0.95 / 0.99（或 90、"95%"、"0.95" 等等价写法）
cs = invert_shift(["1.5", "2", "2", "3.25"], ["2", "2.5"], 0.95)
cs.intervals        # 区间列表，端点为精确 Fraction，None 表示无界
Fraction(0) in cs   # 成员判断：p(theta) 严格大于显著性水平才在集合内
```

- `Interval(lower, upper, lower_closed, upper_closed)`：每个连通分支的开闭边界；
  单点、缺口、开闭端点、无界部分均按原样保留。
- 显著性水平 α = 1 − 置信度；`p == α` 的 θ 不被接受（开边界）。
- 枚举复杂度为 C(n+m, n)，适用于小样本精确计算。

## 运行测试

```sh
python3 -m unittest discover -s tests -v
```

测试含独立的标签枚举参照实现，核对含重复值且组大小不等样本的 U 与精确 p 值；
对置信集合逐一检查跨组差值本身、相邻差值之间及两侧外部的成员关系，
并覆盖 p 恰好等于显著性水平（不接受）与退化单点集合的情形。
