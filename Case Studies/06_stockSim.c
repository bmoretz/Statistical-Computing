void
stockSim( double *ans, const int *len, const double *rho,
          const double *psi, const double *eps)
{
  int i, j;
  
  double psi_rho = (1 - *rho) * (*psi);
  
  for(i = 1, j = *len + 1; i < *len; i++, j++)
  {
    ans[i] = *rho * ans[i - 1] + psi_rho * ans[j - 1] + eps[i];
    ans[j] = *rho * ans[j - 1] + psi_rho * ans[i - 1] + eps[j];
  }
}