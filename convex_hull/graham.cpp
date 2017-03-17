int graham(int n, point * p, point * ch, bool comEdge = false){
    if(n < 3){ for(int i = 0; i < n; i++) ch[i] = p[i]; return n; }
    const double e1 = comEdge ? eps : -eps; int i, j, k;
    sort(p, p + n); ch[0] = p[0]; ch[1] = p[1];
    for(i = j = 2; i < n; ch[j++] = p[i++]) while(j > 1 && (ch[j - 2] - ch[j - 1]) * (p[i] - ch[j - 1]) > e1) j--;
    ch[k = j++] = p[n - 2];
    for(i = n - 3; i > 0; ch[j++] = p[i--]) while(j > k && (ch[j - 2] - ch[j - 1]) * (p[i] - ch[j - 1]) > e1) j--;
    while (j > k && (ch[j - 2] - ch[j - 1]) * (ch[0] - ch[j - 1]) > e1) j--;
    return j;
}//求凸包,p会被打乱顺序,ch为逆时针,comEdge为true时保留共线点,重点会导致不稳定

void graham(PointSet &p, PointSet &ch, bool comEdge = false){
    if (p.size()<3) {
    	for (int i=0;i<p.size();i++)
    		ch[i] = p[i];
    	return; 
    }
    const double e1 = comEdge ? eps : -eps;
    int i, j, k;
    sort(p.begin(), p.end(), );
    ch[0] = p[0]; ch[1] = p[1];
    for (i = j = 2; i < n; ch[j++] = p[i++]) 
    	while(j > 1 && (ch[j - 2] - ch[j - 1]) * (p[i] - ch[j - 1]) > e1)
    		j--;
    ch[k = j++] = p[n - 2];
    for (i = n - 3; i > 0; ch[j++] = p[i--])
    	while(j > k && (ch[j - 2] - ch[j - 1]) * (p[i] - ch[j - 1]) > e1)
    		j--;
    while (j > k && (ch[j - 2] - ch[j - 1]) * (ch[0] - ch[j - 1]) > e1)
    	j--;
}