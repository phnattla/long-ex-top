import math

m_w = 80.385
m_b = 4.18

m_t = 172.5

E_b = -(m_w**2)/(2*m_t) + (m_b**2)/(2*m_t) + 0.5*m_t

print("Predicted E_peak: "+str(E_b)+" GeV")
