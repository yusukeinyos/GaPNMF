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
        static int K;

        static int[] sorted_index;

        static void Main(string[] args)
        {
            double[,] data = new double[100, 2];
            for (int i = 0; i < 70; i++)
            {
                data[i, 0] = 0.5 * i;
                data[i, 1] = bessel_2(data[i, 0], 0.8);
            }
            CsvFileIO.CsvFileIO.WriteData("output.csv", data);
        }

        static void init()
        {
            X = CsvFileIO.CsvFileIO.ReadData("");
            I = X.GetLength(0);
            J = X.GetLength(1);
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
            sorted_index = (int[])Enumerable.Range(0, K - 1);
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
            Array.Sort(Theta_estimated, sorted_index);
            Array.Reverse(Theta_estimated);
            Array.Reverse(sorted_index);    

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
                b = Mt.Factorial(m) * Mt.Gamma(dim + m + 1.0);
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
    }
}
