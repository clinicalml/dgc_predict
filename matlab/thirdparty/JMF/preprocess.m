g2g_ge_sum = sum(g2g_ge, 2);
g2g_ge_idx = find(g2g_ge_sum > 1);
g2g_ppi_sum = sum(g2g_ppi, 2);
g2g_ppi_idx = find(g2g_ppi_sum > 1);
g2g_seq_sum = sum(g2g_seq, 2);
g2g_seq_idx = find(g2g_seq_sum > 1);
g2g_go_sum =  sum(g2g_go, 2);
g2g_go_idx = find(g2g_go_sum > 1);
g2g_idx1 = intersect(g2g_ppi_idx, g2g_ge_idx);
g2g_idx2 = intersect(g2g_seq_idx, g2g_go_idx);
g2g_idx = intersect(g2g_idx1, g2g_idx2);

g2g_seq = g2g_seq(g2g_idx, g2g_idx);
g2g_go = g2g_go(g2g_idx, g2g_idx);
g2g_ppi = g2g_ppi(g2g_idx, g2g_idx);
g2g_ge = g2g_ge(g2g_idx, g2g_idx);

association = association(:, g2g_idx);


d2d(logical(eye(size(d2d))))=0;

g2g_seq(logical(eye(size(g2g_seq))))=0;
g2g_go(logical(eye(size(g2g_go))))=0;
g2g_ppi(logical(eye(size(g2g_ppi))))=0;
g2g_ge(logical(eye(size(g2g_ge))))=0;