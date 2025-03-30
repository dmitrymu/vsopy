# Simplified Differential Photometry

$$\begin{equation}
M = m - k^{'} X - k^{"} X C + T C + Z
\end{equation}$$

$$\begin{equation}
M_t - M_c = (m_t-m_c) - k^{''} X(C_t-C_c) + T(C_t-C_c)
\end{equation}$$

$$\begin{equation}
C = A - B \\
= (a-b) - (k_A^{'} - k_B^{'})X   - (k_A^{''} - k_B^{''}) X C  + (T_A - T_B) C + (Z_A - Z_B)
\end{equation}$$

$$\begin{equation}
\begin{split}
 C_t-C_c = (a_t - b_t) - (a_c - b_c) - (k_A^{''} - k_B^{''}) X (C_t-C_c) + (T_A - T_B)(C_t-C_c) \\
 = c_t - c_c + [(T_A - T_B) - (k_A^{''} - k_B^{''}) X](C_t-C_c) \\
 =  K(c_t - c_c) \\
 = \frac{c_t - c_c}{1-[(T_A - T_B) - (k_A^{''} - k_B^{''}) X]}
\end{split}
\end{equation}$$

$$\begin{equation}
\Delta_c = c_t - c_c \\
\Delta_a = a_t - a_c
\end{equation}$$

$$\begin{equation}
A_t = A_c + \Delta_a  - \frac{k^{''}_{A} X }{K}\Delta_c
+ \frac{T_A}{K}\Delta_c \\
= A_c + \Delta_a  +\frac{T_A - k^{''}_{A} X}{K} \Delta_c
\end{equation}$$

$$\begin{equation}
    C_t = A_t - B_t = (A_c - B_c) + T_{AB} * ((a_t-b_t) - (a_c-b_c)) \\
    A_t = a_t + (A_c-a_c) + T_A * (C_t - (A_c - B_c)) \\
    B_t = A_t - C_t
\end{equation}$$

$$\begin{equation}
  a - b = T_{AB} (A - B) + Z_{AB} \\
  A - a = T_A (A - B) + Z_A \\
  B - b = T_B (A - B) + Z_B
\end{equation}$$
