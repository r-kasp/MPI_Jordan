//#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <sys/resource.h>
#include <sys/time.h>
#include "matrix.h"
#include <unistd.h>
#include <cstring>

#define ok 0

double get_cpu_time()
{
  struct rusage buf;
  getrusage(RUSAGE_THREAD, &buf);
  return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e+6;
}


double get_full_time()
{
  struct timeval buf;
  gettimeofday(&buf, 0);
  return buf.tv_sec + buf.tv_usec / 1.e+6;
}


int main (int argc, char *argv[])
{
    int n, m, p, r, s, k;
    
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size (comm, &p);
    MPI_Comm_rank (comm, &k);
    
    char * filename = nullptr;
    if (!((argc == 5 || argc == 6) 
        && sscanf(argv[1], "%d", &n) == 1 
        && sscanf(argv[2], "%d", &m) == 1 
        && sscanf(argv[3], "%d", &r) == 1 
        && sscanf(argv[4], "%d", &s) == 1))
    {
        printf ("Usage %s n m r s (file) \n", argv[0]);
        MPI_Finalize ();
        return 0;
    }
  
    if (argc == 6)
        filename = argv[5];
       
    if (p >= n / m + (n % m != 0))
      p = n / m + (n % m != 0);
    if (k >= n / m + (n % m != 0))
    {
      MPI_Finalize ();
      return 0;
    }
    
    int blocked_columns_cnt = (n / m) + (n % m != 0);
    int col_cnt = blocked_columns_cnt / p + (k < blocked_columns_cnt % p); //check
          
    double *a = new double[m * blocked_columns_cnt * m * col_cnt];
    if (!a)
    {
        printf("Not enough memory\n");
        MPI_Finalize ();
        return 0;
    }
    //Вектор b
    double *b;
    //Копия вектора b
    double *c;
    
    b = new double[m * blocked_columns_cnt];
    if (!b)
      {
        printf("Not enough memory\n");
        delete [] a;
        MPI_Finalize ();
        return 0;
      }
    c = new double[m * blocked_columns_cnt];
    if (!c)
      {
        printf("Not enough memory\n");
        delete [] a;
        delete [] b;
        MPI_Finalize ();
        return 0;
      }
    double *buf = new double[m * m * blocked_columns_cnt];
    if (!buf)
      {
        printf("Not enough memory\n");
        delete [] a;
        delete [] b;
        delete [] c;
        MPI_Finalize ();
        return 0;
      }
    
    memset (a, 0, m * blocked_columns_cnt * m * col_cnt * sizeof (double));  
    memset (buf, 0, m * m * blocked_columns_cnt * sizeof (double));  
    memset (b, 0, m * blocked_columns_cnt * sizeof (double));  
    memset (c, 0, m * blocked_columns_cnt * sizeof (double));  
    
    if (filename)
    {
      int ret = read_matrix_col (a, n, m, p, k, buf, filename, comm);
      if (ret != SUCCESS)
        {
          switch (ret)
            {
              case ERROR_OPEN:
                  printf("Cannot open %s\n", filename);
                  break;
              case ERROR_READ:
                  printf("Cannot read %s\n", filename);
                  break;
              default:
                  printf("Unknown error %d in file %s\n", ret, filename);
            }
          delete [] a;
          delete [] b;
          delete [] c;
          delete [] buf;
          MPI_Finalize ();
          return 0;
        }
    }
    else
      init_matrix (a, s, n, m, p, k);
      
    init_vector(b, a, n, m, p, k, comm);
    copy_vector(c, b, n);
    print_matrix_col (a, n, m, p, k, buf, r, comm);
    double norm = mpi_norma_matrix (a, m, n, col_cnt, k, p, comm);
    
    double elapsed = MPI_Wtime();
    //SOLVE
    int ret = solve_jordan (a, buf, b, norm, n, m, k, p, comm);
    elapsed = MPI_Wtime() - elapsed;
    
    if (ret != SUCCESS)
      {
        if (k == 0)
          printf ("Can't Solve. The Matrix is Singular\n");
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] buf;
        MPI_Finalize ();
        return 0;
      }
    
    if (filename)
    {
      int ret = read_matrix_col (a, n, m, p, k, buf, filename, comm);
      if (ret != SUCCESS)
        {
          switch (ret)
            {
              case ERROR_OPEN:
                  printf("Cannot open %s\n", filename);
                  break;
              case ERROR_READ:
                  printf("Cannot read %s\n", filename);
                  break;
              default:
                  printf("Unknown error %d in file %s\n", ret, filename);
            }
          delete [] a;
          delete [] b;
          delete [] c;
          delete [] buf;
          MPI_Finalize ();
          return 0;
        }
    }
    else
      init_matrix (a, s, n, m, p, k);
    
    double residual = count_residual (n, m, a, b, c, k, p, comm); 
    
    if (k == 0)
      {
        printf ("Result : \n");
        print_vector (b, n, r);
        printf ("%s : residual = %e elapsed = %.2f s = %d n = %d m = %d p = %d\n", argv[0], residual, elapsed, s, n, m, p);
      }
    
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] buf;
    MPI_Finalize ();
    return 0;
}
