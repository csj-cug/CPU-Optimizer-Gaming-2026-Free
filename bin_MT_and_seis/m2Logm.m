function ml = m2Logm(m,m_min,m_max)
    ml=log(m-m_min)-log(m_max-m);