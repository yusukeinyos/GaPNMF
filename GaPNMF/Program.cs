using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InoueLab;
using CsvFileIO;

namespace GaPNMF
{
    class Program
    {
        static double[,] X;
        static double[,] W_estimated;
        static double[,] H_estimated;
        static double[] Theta_estimated;

        static double[,] lo_W;
        static double[,] tau_W;
        static double[,] lo_H;
        static double[,] tau_H;
        static double[] lo_Theta;
        static double[] tau_Theta;

        static double[, ,] fai;
        static double[,] omega;

        static double a;
        static double b;
        static double c;
        static double alpha;

        static int I;
        static int J;
        static int K; //最初に用意しておく基底
        static int L; //実際に使う基底

        static int[] goodness_index;

        static void Main(string[] args)
        {
            //double[,] data = new double[100, 2];
            //for (int i = 0; i < 70; i++)
            //{
            //    data[i, 0] = 0.5 * i;
            //    data[i, 1] = bessel_2(data[i, 0], 0.8);
            //}
            //CsvFileIO.CsvFileIO.WriteData("output.csv", data);

            int max_itteration = 50;
            int itteration = 0;
            init();
            double[,] theta_regist = new double[max_itteration, K];
            do
            {
                Update();
                Console.WriteLine("itteration : " + itteration + 1);
                for (int k = 0; k < K; k++)
                    theta_regist[itteration, k] = GIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]);
                itteration++;
            } while (itteration < max_itteration);
            CsvFileIO.CsvFileIO.WriteData("theta_regist.csv", theta_regist);

            //double[,] data = new double[100,1];
            //for (int i = 0; i < 100; i++)
            //    data[i,0] = gamma_rnd(2, 1);
            //CsvFileIO.CsvFileIO.WriteData("test.csv", data);

        }

        static void init()
        {
            X = CsvFileIO.CsvFileIO.ReadData(@"C:\Users\優\Desktop\音素材\mix4.csv");
            I = X.GetLength(0);
            J = X.GetLength(1);
            K = 10;
            L = K;
            lo_W = new double[I, K];
            tau_W = new double[I, K];
            W_estimated = new double[I, K];
            lo_H = new double[K, J];
            tau_H = new double[K, J];
            H_estimated = new double[K, J];
            lo_Theta = new double[K];
            tau_Theta = new double[K];
            Theta_estimated = new double[K];
            fai = new double[I, J, K];
            omega = new double[I, J];

            for (int k = 0; k < K; k++)
            {
                for (int i = 0; i < I; i++)
                    tau_W[i, k] = 0.1;
                for (int j = 0; j < J; j++)
                    tau_H[k, j] = 0.1;
                tau_Theta[k] = 0.1;
                lo_Theta[k] = initMatrix(1, 1)[0, 0];
            }

            lo_W = initMatrix(I, K);
            lo_H = initMatrix(K, J);


            a = 0.1;
            b = 0.1;
            c = 1.0;
            alpha = 1.0;
        }

        static void estimated()
        {
            List<int> threshold = new List<int>();
            for (int k = 0; k < K; k++)
                if (GIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]) > 0.01)
                    threshold.Add(k);
            int L = threshold.Count;
            W_estimated = new double[I, L];
            H_estimated = new double[L, J];

        }

        static void Update()
        {
            updateFai();
            updateOmega();
            updateHpara();
            updateWpara();
            updateThetapara();
        }

        static void updateWpara()
        {
            for (int i = 0; i < I; i++)
                for (int k = 0; k < K; k++)
                {
                    double sum1 = 0, sum2 = 0;
                    for (int j = 0; j < J; j++)
                    {
                        sum1 += GIG_expectation(b, lo_H[k, j], tau_H[k, j]) / omega[i, j];
                        sum2 += X[i, j] * fai[i, j, k] * fai[i, j, k] * invGIG_expectation(b, lo_H[k, j], tau_H[k, j]);
                    }

                    lo_W[i, k] = a + GIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]) * sum1;
                    tau_W[i, k] = invGIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]) * sum2;
                }
        }

        static void updateHpara()
        {
            for (int k = 0; k < K; k++)
                for (int j = 0; j < J; j++)
                {
                    double sum1 = 0, sum2 = 0;
                    for (int i = 0; i < I; i++)
                    {
                        sum1 += GIG_expectation(a, lo_W[i, k], tau_W[i, k]) / omega[i, j];
                        sum2 += X[i, j] * fai[i, j, k] * fai[i, j, k] * invGIG_expectation(a, lo_W[i, k], tau_W[i, k]);
                    }

                    lo_H[k, j] = b + GIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]) * sum1;
                    tau_H[k, j] = invGIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]) * sum2;
                }
        }

        static void updateThetapara()
        {
            for (int k = 0; k < K; k++)
            {
                double sum1 = 0, sum2 = 0;
                for (int i = 0; i < I; i++)
                    for (int j = 0; j < J; j++)
                    {
                        sum1 += GIG_expectation(a, lo_W[i, k], tau_W[i, k]) * GIG_expectation(b, lo_H[k, j], tau_H[k, j]) / omega[i, j];
                        sum2 += X[i, j] * fai[i, j, k] * fai[i, j, k] * invGIG_expectation(a, lo_W[i, k], tau_W[i, k]) * invGIG_expectation(b, lo_H[k, j], tau_H[k, j]);
                    }

                lo_Theta[k] = alpha * c + sum1;
                tau_Theta[k] = sum2;
            }
        }

        static void updateFai()
        {
            for (int k = 0; k < K; k++)
            {
                double sum = 0;
                for (int i = 0; i < I; i++)
                    for (int j = 0; j < J; j++)
                    {
                        fai[i, j, k] = 1.0 / (invGIG_expectation(a, lo_W[i, k], tau_W[i, k]) * invGIG_expectation(b, lo_H[k, j], tau_H[k, j]) * invGIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]));
                        sum += fai[i, j, k];
                    }
                for (int i = 0; i < I; i++)
                    for (int j = 0; j < J; j++)
                        fai[i, j, k] /= sum;        //kについて規格化
            }
        }

        static void updateOmega()
        {
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < K; k++)
                        sum += GIG_expectation(a, lo_W[i, k], tau_W[i, k]) * GIG_expectation(b, lo_H[k, j], tau_H[k, j]) * GIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]);
                    omega[i, j] = sum;
                }
        }

        static void sortingOnTheta()
        {
            //降順にソート
            for (int k = 0; k < K; k++)
                Theta_estimated[k] = GIG_expectation(alpha / K, lo_Theta[k], tau_Theta[k]);

            goodness_index = (int[])Theta_estimated.Select((num, index) => new { Index = index, Value = num }).Where(num => num.Value > 0.001).Select(num => num.Index);
            L = Theta_estimated.Where(num => num > 0.001).Count();


        }

        //--------------------------------------------------------------------
        static double GIG_expectation(double gamma, double lo, double tau)
        {
            double a = Math.Sqrt(lo * tau);
            double b = bessel_2(2 * a, gamma + 1.0);
            double c = bessel_2(2 * a, gamma);
            return b * Math.Sqrt(tau) / (c * Math.Sqrt(lo));

        }

        static double invGIG_expectation(double gamma, double lo, double tau)
        {
            double a = Math.Sqrt(lo * tau);
            double b = bessel_2(2 * a, gamma - 1.0);
            double c = bessel_2(2 * a, gamma);
            return b * Math.Sqrt(lo) / (c * Math.Sqrt(tau));
        }

        static double bessel_1(double x, double dim)
        {
            double y = 0;
            double a = Math.Pow(x / 2.0, dim);
            double b = 0;
            double delta_y = 0;
            int m = 0;
            do
            {
                b = Mt.Factorial(m) * Gamma2(dim + m + 1.0);
                delta_y = 1.0 / b * Math.Pow(x / 2.0, 2.0 * m);
                if (m % 2 == 1)
                    delta_y = -delta_y;

                y += delta_y;
                m++;
            } while (Math.Abs(delta_y) > 0.000000000000000001);
            y = y * a;
            return y;
        }

        static double bessel_2(double x, double dim)
        {
            double y = 0;
            double a = bessel_1(x, dim);
            double b = bessel_1(x, -dim);
            y = (a * Math.Cos(dim * Math.PI) - b) / Math.Sin(dim * Math.PI);
            return y;
        }

        static double Gamma2(double d)
        {
            if (d > 0)
                return Mt.Gamma(d);
            else
                return Math.PI / (Math.Sin(Math.PI * d) * Mt.Gamma(1.0 - d));

        }
        static double gamma_rnd(double shape, double scale)
        {
            double output;
            RandomMT rand = new RandomMT();
            double input = rand.Double();
            output = scale * Mt.InverseIncompleteGamma(shape, input);
            return output;
        }

        //初期値をランダムに設定 (0,1]
        public static double[,] initMatrix(int N, int M)
        {
            RandomMT rand = new RandomMT();
            double[,] A = new double[N, M];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    A[i, j] = 10 * rand.Double32OC();
            return A;
        }

        static double[,] pivotMatrix(double[,] mat, int[] index_key, int flag)
        {
            int I = mat.GetLength(0);
            int J = mat.GetLength(1);
            double[] re;
            if (flag == 0)          //行ピボット
            {
                re = new double[J];
                for (int i = 0; i < I; i++)
                    for (int j = 0; j < J; j++)
                    {
                        re[j] = mat[i, j];
                        mat[i, j] = mat[index_key[i], j];
                        mat[index_key[i], j] = re[j];
                    }
            }
            else if (flag == 1)     //列ピボット
            {
                re = new double[I];
                for (int j = 0; j < J; j++)
                    for (int i = 0; i < I; i++)
                    {
                        re[i] = mat[i, j];
                        mat[i, j] = mat[i, index_key[j]];
                        mat[i, index_key[j]] = re[i];
                    }
            }
            return mat;
        }
    }
}
