//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include "matrix.h"

static int loc2glob (int m, int loc, int p, int k)
{
  return (loc / m) * m * p + k * m + (loc % m);
}

static void fill_arranged_items(int *array, int size)
{
	for (int i = 0; i < size; i++)
		array[i] = i;
}


static void set_block(double *a, int i, int j, double *b, int n, int m, int height, int width)
{
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			 a[(i*m+h)*n + j*m + w] = b[h*width + w];
		}
	}
}



static void get_block (double *a, int i, int j, double *b, int n, int m, int height, int width)
{
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
		  b[h * width + w] = a[(i * m + h) * n + j * m + w];
		}
	}
}

 
static void get_vector(double *b, double *res, int ind, int m, int size)

{
	for (int i = 0; i < size; i++)
		res[i] = b[ind*m+i];
}


static void set_vector(double *b, double *buf, int ind, int m, int size)
{
	for (int i = 0; i < size; i++)
		b[ind*m + i] = buf[i];
}


void subtract_blocks(double *a, double *b, int n, int m)
{
	for (int i = 0; i < n*m; i++)
		a[i] -= b[i];
}


void clean_memory(double *buf, double *buf_mm, double *buf_mm2, double *buf_ostost, double *buf_ostost2, double *buf_ost, double *buf_m, double *buf_m2, int *r_m, int *r_ost)
{
  delete [] buf;
	delete [] buf_mm;
	delete [] buf_mm2;
	delete [] buf_ostost;
  delete [] buf_ostost2;
	delete [] buf_ost;
	delete [] buf_m;
	delete [] buf_m2;
	delete [] r_m;
	delete [] r_ost;
}


static void multiply_block_on_str (int i, int start, double *a, double *mult, int col_cnt, int m, double *buf, double *res)
{
  for (int w = start; w < col_cnt; w++)
    {
	    get_block (a, i, w, buf, col_cnt * m, m, m, m);
	    multiply_blocks_m (mult, buf, res, m);
	    set_block (a, i, w, res, col_cnt * m, m, m, m);
    }
}

static void multiply_block_on_vector (double *buf, double *b, int str_ind, int m, int block_sz, double *b_old)
{
	copy_vector(b_old, b + str_ind * m, block_sz);
	for (int i = 0; i < block_sz; i++)
	  {
		  double s = 0;
		  for (int j = 0; j < block_sz; j++)
		    {
			    s += buf[i * block_sz + j] * b_old[j];
		    }
		  b[str_ind * m + i] = s;
	  }
}

static void subtract_with_multiply (double *a, int from_ind, int ind, int start_stl, double *mult, int col_cnt, int m, double *buf, double *res)
{
	for (int w = start_stl; w < col_cnt; w++)
		{
			get_block (a, ind, w, buf, col_cnt * m, m, m, m);
			multiply_blocks_m (mult, buf, res, m);
			get_block (a, from_ind, w, buf, col_cnt * m, m, m, m);
			subtract_blocks (buf, res, m, m);
			set_block (a, from_ind, w, buf, col_cnt * m, m, m, m);
		}
}

static void change_str_blocks(int start, double *a, double *b, int k, int ind_max, int m, int col_cnt, double *buf_mm, double *buf_mm2, double *buf_m, double *buf_m2)
{
	for (int i = start; i < col_cnt; i++)
	{
	  get_block (a, k, i, buf_mm, col_cnt * m, m, m, m);
	  get_block (a, ind_max, i, buf_mm2, col_cnt * m, m, m, m);
	  
	  set_block (a, k, i, buf_mm2, col_cnt * m, m, m, m);
	  set_block (a, ind_max, i, buf_mm, col_cnt * m, m, m, m);
	}
  get_vector (b, buf_m, k, m, m);
  get_vector (b, buf_m2, ind_max, m, m);
  
  set_vector (b, buf_m2, k, m, m);
  set_vector (b, buf_m, ind_max, m, m);
}

static void subtract_with_multiply_vector (double *b, int from_ind, int ind, double *mult, int m, double *buf, double *res)
{
  get_vector (b, buf, ind, m, m);
	multiply_matrix_and_vector (mult, buf, res, m, m);
	get_vector (b, buf, from_ind, m, m);
	subtract_blocks (buf, res, m, 1);
	set_vector (b, buf, from_ind, m, m); 
}

static void set_unitary_matrix (double *a, int n)
{
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        {
          a[i * n + j] = (i == j);
        }
    }
}

static int get_reversed_matrix_for_block(double *a, int n, int m, int size, int i, int j, double *res, double eps, double *buf, int * r)
{
	get_block (a, i, j, buf, n, m, size, size);
	fill_arranged_items (r, size);
	set_unitary_matrix (res, size);
	int ret = get_reversed_matrix_gauss (size, buf, res, r, eps); 
	return ret;
}

static int get_main_block (double *a, int n, int m, int k, int p, int col_cnt, int stl, double eps, double *buf, double *buf2, int *r)
{
  int s = n / m;
  double min_norm = 1e308;
  int res = -1;
  (void)k;
  for (int i = stl; i < s; i++)
  {
    int ret = get_reversed_matrix_for_block (a, col_cnt * m, m, m, i, stl / p, buf, eps, buf2, r);
    if (ret == SING_MATRIX)
      continue;
    double check = norma_matrix (buf, m, m);
    if (check < min_norm)
    {
      min_norm = check;
      res = i;
    }
  }
  return res;
}

static void copy_column (int column, int col_cnt, double *a_buf, double *a, int m, int height)
{
  for (int i = 0; i < height; i++)
    {
      for (int j = 0; j < m; j++)
        {
          a_buf[i * m + j] = a[i * col_cnt * m + column * m + j];
        }
    }
}

int solve_jordan (double *a, double *a_buf, double *b, double norm, int n, int m, int k, int p, MPI_Comm comm)
{
  int s = n / m;
	int ost = n % m;
	double eps = 1e-16 * norm;// 1e-14
	int blocked_columns_cnt = (n / m) + (n % m != 0);
  int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p);
	
	double *buf = new double[m*m];
  double *buf_mm = new double[m*m];
	double *buf_mm2 = new double[m*m];
	double *buf_ostost = new double[ost*ost];
  double *buf_ostost2 = new double[ost*ost];
	double *buf_ost = new double[ost];
	double *buf_m = new double[m];
	double *buf_m2 = new double[m];
	int *r_m = new int[m];
	int *r_ost = new int[ost];
	
	for (int step = 0; step < s + (ost != 0); step++)
	  {
	    int owner = step % p;
	    int ret = SUCCESS;
	    
	    //Получаем текущий стобец
	    if (owner == k)
	      copy_column (step / p, col_cnt, a_buf, a, m, blocked_columns_cnt * m);
	    MPI_Bcast (a_buf, m * m * blocked_columns_cnt, MPI_DOUBLE, owner, comm);
	    
	    if (owner == k)
	      {
	        if (ost != 0 && step == s) //если последний столбец, который не поделился
	          {
	            get_block (a, s, step / p, buf, col_cnt * m, m, ost, ost);
	            if (norma_matrix (buf, ost, ost) < eps)
	              ret = SING_MATRIX;
	            else
	              ret = s;
	          }
	        else
	          ret = get_main_block (a, n, m, k, p, col_cnt, step, eps, buf, buf_mm, r_m);//поиск главного элемента
	      }
	      
	    MPI_Bcast (&ret, 1, MPI_INT, owner, comm);
	    
	    if (ret < 0)
	      {
	        clean_memory (buf, buf_mm, buf_mm2, buf_ostost, buf_ostost2, buf_ost, buf_m, buf_m2, r_m, r_ost);
	        return ret;
	      }
	    
	    if (owner == k)
	      {
	        if (ost != 0 && step == s)
	          {
	            get_reversed_matrix_for_block (a, col_cnt * m, m, ost, ret, s / p, buf_ostost, eps, buf_ostost2, r_ost);
	            memset (buf, 0, m * m  * sizeof (double)); 
	            set_block (buf, 0, 0, buf_ostost, m, m, ost, ost);
	          }
	        else
	          get_reversed_matrix_for_block (a, col_cnt * m, m, m, ret, step / p, buf, eps, buf_mm, r_m);
	      }
	    
	    MPI_Bcast (buf, m * m, MPI_DOUBLE, owner, comm);
	    
	    //Получили главный элемент
	    multiply_block_on_str (ret, step / p + (owner == k), a, buf, col_cnt, m, buf_mm, buf_mm2);
	    multiply_block_on_vector (buf, b, ret, m, m, buf_m);
	  
	    for (int i = 0; i < s + (ost != 0); i++)
	      {
	        if (i == ret)
	          continue;
	        get_block (a_buf, i, 0, buf, m, m, m, m);
	        if (norma_matrix (buf, m, m) < eps)
	          continue;
	        subtract_with_multiply (a, i, ret, step / p + (owner == k), buf, col_cnt, m, buf_mm, buf_mm2); //CHECK
	        subtract_with_multiply_vector (b, i, ret, buf, m, buf_m, buf_m2);
	      }
	    if (ret != step)
	      change_str_blocks (step / p + (owner == k), a, b, step, ret, m, col_cnt, buf_mm, buf_mm2, buf_m, buf_m2);	   
	  }
    
  clean_memory (buf, buf_mm, buf_mm2, buf_ostost, buf_ostost2, buf_ost, buf_m, buf_m2, r_m, r_ost);
  return SUCCESS;
}

void mpi_multiply_matrix_on_vector (double *a, double *b, double *res, int n, int m, int k, int p, MPI_Comm comm)
{
  int blocked_columns_cnt = (n / m) + (n % m != 0);
  int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p);
  double var = 0;
  for (int i = 0; i < n; i++)
    {
      double s = 0;
      for (int j_loc = 0; j_loc < col_cnt * m; j_loc++)
        {
          int j_glob = loc2glob (m, j_loc, p, k);
          if (j_glob < n)
            s += a[i * col_cnt * m + j_loc] * b[j_glob];
        }
      res[i] = s;
      var = s;
      for (int process = 0; process < p; process++)
        if (process != k)
          MPI_Send (&var, 1, MPI_DOUBLE, process, 0, comm);
      MPI_Status st;
      for (int process = 0; process < p; process++)
        {
          if (process != k)
            {
              MPI_Recv (&var, 1, MPI_DOUBLE, process, 0, comm, &st);
              res[i] += var;
            }
        }     
    }
}
 

double count_residual (int n, int m, double *a, double *b, double *c, int k, int p, MPI_Comm comm)
{
  double *rsd = new double[n];
  mpi_multiply_matrix_on_vector (a, b, rsd, n, m, k, p, comm);
  subtract_vectors (rsd, c, n);
  double norm = norma_vector (c, n);
  if (norm >= 0 && norm <= 0)
    {
      delete [] rsd;
      return -1;
    }
  double res = norma_vector (rsd, n) / norm;
  delete [] rsd;
  return res;
}


static int read_array (FILE * fp, double * a, int len)
{
  for (int i = 0; i < len; i++)
  {
    if (fscanf (fp, "%lf", a + i) != 1)
      return ERROR_READ;
  }
  return SUCCESS;
}


int read_matrix_col(double *a, int n, int m, int p, int k, double *buf, char *name, MPI_Comm comm)
{
  int main_k = 0;
  int err = SUCCESS;
  FILE *fp = nullptr;
  
  if (k == main_k)
    {
      fp = fopen (name, "r");
      if (fp == nullptr)
        err = ERROR_OPEN;
    }
    
  MPI_Bcast (&err, 1, MPI_INT, main_k, comm);
  
  if (err)
    return err;
    
  double *block_buf = new double[m * m];
  if (!block_buf)
    return ERROR_MEMORY;
    
  int b, max_b = (n + m - 1) / m;
  
  int row_size = (m <= n - (max_b - 1) * m ? m : n - (max_b - 1) * m);
  memset (buf, 0, n * row_size * sizeof (double));
  
  int blocked_columns_cnt = (n / m) + (n % m != 0);
  int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p);
  
  for (b = 0; b < max_b; b++)
    {
      int rows = (m <= n - b * m ? m : n - b * m);
      int height = m;
      if (n % m != 0 && b == max_b - 1)
        height = n % m;
      if (k == main_k)
        {
          err += read_array (fp, buf, n * rows);
          for (int j = 0; j < max_b; j++)
            {
              int owner = j % p;
              int width = m;
              int j_loc = j / p;
              if (n % m != 0 && j == max_b - 1)
                  width = n % m;
              get_block (buf, 0, j, block_buf, n, m, height, width);
              
              if (owner == k)
                {
                  set_block (a, b, j_loc, block_buf, col_cnt * m, m, height, width);
                }
              else
                {
                  MPI_Send (block_buf, height * width, MPI_DOUBLE, owner, 0, comm);
                }
            }
        }
      else
        {
          for (int j = 0; j < max_b; j++)
            {
              int owner = j % p;
              int width = m;
              int j_loc = j / p;
              if (n % m != 0 && j == max_b - 1)
                  width = n % m;
              if (owner == k)
              {
                MPI_Status st;
                MPI_Recv (block_buf, height * width, MPI_DOUBLE, main_k, 0, comm, &st);
                set_block (a, b, j_loc, block_buf, col_cnt * m, m, height, width);
              }
            }
        }
    }
  if (k == main_k)
  {
    fclose (fp);
    fp = nullptr;
  }
  MPI_Bcast (&err, 1, MPI_INT, main_k, comm);
  delete [] block_buf;
  if (err != SUCCESS)
    return err;
  return SUCCESS;
}



double f(int s, int n, int i, int j)
{
    i++; j++;
    if (s == 1)
        return n - (i > j ? i : j) + 1;
    if (s == 2)
        return (i > j ? i : j);
    if (s == 3)
        return (i > j ? i-j : j-i);
    if (s == 4)
        return (double)1/(i+j-1); 
    return 1;
}

void init_matrix(double *a, int s, int n, int m, int p, int k)
{
  int blocked_columns_cnt = (n / m) + (n % m != 0);
  int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p);
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < col_cnt * m; j++)
        {
          int global_j = (j / m) * m * p + k * m + (j % m);
          a[i * col_cnt * m + j] = f (s, n, i, global_j);
        }
    }
  
}


void print_matrix_col (double *a, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm comm)
{
  int main_k = 0;
  int p_n = (max_print < n ? max_print : n);
  int max_b = (n + m - 1) / m;
  int blocked_columns_cnt = (n / m) + (n % m != 0);
  int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p);
  
  for (int i = 0; i < p_n; i++)
    {
      for (int j = 0; j < max_b; j++)
        {
          int len = ((n % m != 0 && j == max_b - 1) ? n % m : m);
          int owner = j % p;
          int j_loc = j / p;
          if (k == main_k)
            {
              if (owner == main_k)
                {
                  for (int q = 0; q < len && j * m + q < max_print; q++)
                    printf (" %10.3e", (a + i * m * col_cnt + j_loc * m)[q]);
                }
              else
                {
                  MPI_Status st;
                  MPI_Recv (buf, len, MPI_DOUBLE, owner, 0, comm, &st);
                  for (int q = 0; q < len && j * m + q < max_print; q++)
                    printf (" %10.3e", buf[q]);
                }
            }
          else
            {
              if (owner == k)
                {
                  MPI_Send (a + i * m * col_cnt + j_loc * m, len, MPI_DOUBLE, main_k, 0, comm);
                }
            }
        }
      if (k == main_k)
        printf ("\n");
    }
}


double norma_vector (double *b, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
    	double check = fabs(b[i]);
        if (check > res)
        	res = check;
    }
    return res;
}


double mpi_norma_matrix (double *a, int m, int n, int col_cnt, int k, int p, MPI_Comm comm)
{
  double res = 0;
  for (int i = 0; i < m; i++)
  {
    double check = 0;
    for (int j = 0; j < col_cnt * m; j++)
    {
      if (loc2glob (m, j, p, k) < n)
        check += fabs (a[i * col_cnt * m + j]);
      else
        break;
    }
    for (int process = 0; process < p; process++)
      if (process != k)
        MPI_Send (&check, 1, MPI_DOUBLE, process, 0, comm);
    MPI_Status st;
    double var = 0;
    for (int process = 0; process < p; process++)
      {
        if (process != k)
          {
            MPI_Recv (&var, 1, MPI_DOUBLE, process, 0, comm, &st);
            check += var;
          }
      }     
    if (check > res)
    {
      res = check;
    }
  }
  return res;
}

//horizontal norm
double norma_matrix (double *a, int m, int n)
{
  double res = 0;
  for (int i = 0; i < m; i++)
  {
    double check = 0;
    for (int j = 0; j < n; j++)
    {
      check += fabs(a[i*n + j]);
    }
    if (check > res)
    {
      res = check;
    }
  }
  return res;
}


void multiply_matrix(double *a, double *b, double *c, int m, int n, int k)
{
    for (int i = 0; i < m; i++)
    {
        double * cbuf = c + i * k;
        for (int j = 0; j < k; j++)
          cbuf[j] = 0;
        for (int j = 0; j < n; j++)
        {
          double mult = a[i * n + j];
          double *bbuf = b + j * k;
          for (int q = 0; q < k; q++)
            cbuf[q] += mult * bbuf[q];
        }
    }
}


void multiply_blocks_m (double *a, double *b, double *c, int m)
{
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    int to = (m / 3) * 3;
    //идём до последней строчки которая не поделилась
    for (int i = 0; i < to; i += 3)
    {
      //сначала по обычным блокам
      for (int j = 0; j < to; j += 3)
      {
        s00 = 0; 
        s01 = 0; 
        s02 = 0;
        s10 = 0; 
        s11 = 0; 
        s12 = 0;
        s20 = 0; 
        s21 = 0; 
        s22 = 0;
        for (int q = 0; q < m; q++)
        {
          double a00 = a[i * m + q], a10 = a[(i + 1) * m + q], a20 = a[(i + 2) * m + q],
                 b00 = b[q * m + j], b01 = b[q * m + j + 1], b02 = b[q * m + j + 2];
          s00 += a00 * b00; 
          s01 += a00 * b01; 
          s02 += a00 * b02;
          s10 += a10 * b00; 
          s11 += a10 * b01; 
          s12 += a10 * b02;
          s20 += a20 * b00; 
          s21 += a20 * b01; 
          s22 += a20 * b02;
        }
        c[i * m + j] = s00; 
        c[i * m + j + 1] = s01; 
        c[i * m + j + 2] = s02;
        c[(i + 1) * m + j] = s10; 
        c[(i + 1) * m + j + 1] = s11; 
        c[(i + 1) * m + j + 2] = s12;
        c[(i + 2) * m + j] = s20; 
        c[(i + 2) * m + j + 1] = s21; 
        c[(i + 2) * m + j + 2] = s22;
      }
      //теперь остался один неполный блок
      for (int j = to; j < m; j++)
      {
        s00 = 0; 
        s10 = 0; 
        s20 = 0;
        for (int q = 0; q < m; q++)
        {
          double b00 = b[q * m + j];
          double a00 = a[i * m + q], a10 = a[(i + 1) * m + q], a20 = a[(i + 2) * m + q];
          s00 += a00 * b00;
          s10 += a10 * b00;
          s20 += a20 * b00;
        }
        c[i * m + j] = s00;
        c[(i + 1) * m + j] = s10;
        c[(i + 2) * m + j] = s20;
      }
    } 
    //последняя не поделившаяся строчка
    for (int i = to; i < m; i++)
    {
      //идем до маленького блока
      for (int j = 0; j < to; j += 3)
      {
        s00 = 0; 
        s01 = 0; 
        s02 = 0;
        for (int q = 0; q < m; q++)
        {
          double a00 = a[i * m + q];
          double b00 = b[q * m + j], b01 = b[q * m + j + 1], b02 = b[q * m + j + 2];
          s00 += b00 * a00;
          s01 += b01 * a00;
          s02 += b02 * a00; 
        }
        c[i * m + j] = s00;
        c[i * m + j + 1] = s01;
        c[i * m + j + 2] = s02;
      }
      //маленький блок
      for (int j = to; j < m; j++)
      {
        s00 = 0;
        for (int q = 0; q < m; q++)
          s00 += a[i * m + q] * b[q * m + j];
        c[i * m + j] = s00;
      }
    }
}


void multiply_matrix_and_vector(double *a, double *b, double *c, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
    	  double s = 0;
        for (int j = 0; j < n; j++)
        {
            s += a[i*n + j] * b[j];
        }
        c[i] = s;
    }
}

void subtract_vectors(double *b, double *c, int n)
{
	for (int i = 0; i < n; i++)
		b[i] -= c[i];
}

void copy_matrix(double *x, double *a, int n)
{
	for (int i = 0; i < n*n; i++)
			x[i] = a[i];
}

void copy_vector(double *c, double *b, int n)
{
	for (int i = 0; i < n; i++)
		c[i] = b[i];
}

void init_vector(double *b, double *a, int n, int m, int p, int k, MPI_Comm comm)
{
  int blocked_columns_cnt = (n / m) + (n % m != 0);
  int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p);
  double *var = new double[1];
  
	for (int i = 0; i < n; i++)
	  {
		  double s = 0;
		  for (int j_loc = 0; j_loc < col_cnt * m; j_loc++)
		    {
		      int j_glob = loc2glob (m, j_loc, p, k);
		      if (j_glob % 2 == 0 && j_glob < n)
		        s += a[i * col_cnt * m + j_loc];
		    }
      b[i] = s;
      var[0] = s;
      for (int process = 0; process < p; process++)
        if (process != k)
          MPI_Send (var, 1, MPI_DOUBLE, process, 0, comm);
      MPI_Status st;
      for (int process = 0; process < p; process++)
        {
          if (process != k)
            {
              MPI_Recv (var, 1, MPI_DOUBLE, process, 0, comm, &st);
              b[i] += var[0];
            }
        }
	  }
	  delete [] var;
}

void print_vector(double *b, int n, int r)
{
	double n_max;
  n_max = (n > r ? r : n);
  for (int i = 0; i < n_max; i++)
  	printf(" %10.3e", b[i]);
  printf("\n");
}

//меняет столбцы
static void changestl(double *a, int n, int i, int j)
{
	int q; double c;
	for (q = 0; q < n; q++)
	{
		c = a[q*n+i]; a[q*n+i] = a[q*n+j]; a[q*n+j] = c;
	}
}

//меняет строки начиная с номера start в каждой строке
static void changestr(double *a, int n, int i, int j, int start)
{
	int q; double c; int in = i*n, jn = j*n;
	for (q = start; q < n; q++)
	{
		c = a[in+q]; a[in+q] = a[jn+q]; a[jn+q] = c;
	}
}

//отнимает от i-ой строки j-ую, умноженную на d, начиная с элемента start
static void operation(double *a, int n, int i, int j, double d, int start)
{
	int q; int in = i*n, jn = j*n;
	for (q = start; q < n; q++)
		a[in+q] -= d*a[jn+q];
}


//эта функция с 1 курса, поэтому написана не очень красиво
int get_reversed_matrix_gauss(int n, double *a, double *x, int *r, double eps)
{
	int i,j,imax = 0,jmax = 0,c,kn,in;
	double MAX = 0, aij;
	int k = 0;
	//прямой ход
	while (k < n)
	{
		for (i = k; i < n; i++)
		{
			in = i*n;
			for (j = k; j < n; j++)
			{
				aij = fabs(a[in+j]);
				if (aij > MAX)
				{
					MAX = aij; imax = i; jmax = j;
				}
			}
		}
		if (MAX <= eps)
			return SING_MATRIX;
		if (imax != k)
		{
			changestr(a,n,k,imax,k);
			changestr(x,n,k,imax,0);
		}
		if (jmax != k)
		{
			i = r[k];
			r[k] = r[jmax];
			r [jmax] = i;
			changestl(a,n,k,jmax);
		}
		kn = k*n;
		//проверка что aij не равен нулю была выше
		aij = a[kn+k];
		//делим строки
		for (i = k+1; i < n; i++)
			a[kn+i] /= aij;
		for (i = 0; i < n; i++)
			x[kn+i] /= aij;
		for (i = k+1; i < n; i++)
		{
			aij = a[i*n+k];
			operation(a,n,i,k,aij,k+1);
			operation(x,n,i,k,aij,0);
		}
		MAX = 0;
		k++;
	}
	k = n-1;
	//обратный ход
	while (k > 0)
	{
		for (i = k-1; i >= 0; i--)
		{
			aij = a[i*n+k];
			operation(x,n,i,k,aij,0);
		}
		k--;
	}
	//меняем "столбцы" обратно
	//меняет за квадрат
	
	for (i = 0; i < n; i++)
	{
		c = r[i];
		while (c != i)
		{
			changestr(x,n,c,i,0);
			j = r[c];
			r[c] = r[i];
			r[i] = j;
			c = r[i];
		}
	}
	return SUCCESS;
}
