vip_null_R <- function(x, y, perm, lv_for_vip){
  m <- nrow(y); n <- ncol(y);
	tmp1 <- matrix(0, ncol(x), perm)
	library(pracma)
	
	for (i in 1:perm){
		ind <- sample(m, m, replace=FALSE)
		X1_r = x[ind,]
		result_pls_rand = pls_R(X1_r, y,lv_for_vip)
		
		w <- result_pls_rand$pls_w[,1:lv_for_vip]
		r2 <- result_pls_rand$pls_ssq[1:lv_for_vip,4]
		result_vip_rand = vip_R_v2(w, r2)
		tmp1[,i] = result_vip_rand
		cat(paste0(i, ' '))
				
	}
	tmp1 <- Reshape(tmp1, ncol(x)*perm, 1)
	return(tmp1)

}