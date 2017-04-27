package matrprogonka;

import java.util.Scanner;

public class MatrProgonka {

    static public void logic() {
        Scanner sc = new Scanner(System.in);
        System.out.println("Введите m размерность");
        int m;
        m = sc.nextInt();
        System.out.println("Введите a по оси X");
        double a = sc.nextDouble();
        System.out.println("Введите M количество слоев");
        double M = sc.nextDouble();
        System.out.println("Введите T по оси Y");
        double T = sc.nextDouble();
        //  System.out.println("Введите eps");
        double eps;
        eps = 0.0001;
        //eps = sc.nextDouble();
        double tau = (double) (T / M);
        double z;
        double[] lam = new double[3];
        double[] mu = new double[3];
        double t;
        double gamma;
        double[] x = new double[m + 1];
        double[] A = new double[m + 1];
        double[] B = new double[m + 1];
        double[] C = new double[m + 1];
        double[] F = new double[m + 1];
        double[] u_x = new double[m + 1];
        double[] alfa = new double[m + 1];
        double[] beta = new double[m + 1];
        double[] u = new double[m + 1];
        double[] Yn = new double[m + 1];
        double h = a / m;

        for (int j = 0; j < M; j++) {
            t = tau * j;
            System.out.println("Слой: " + j);
            //null sloi
            /*-----------------*/
            for (int i = 0; i < m + 1; i++) {
                x[i] = i * h;
                u_x[i] = Math.pow(x[i], 3) + x[i] * Math.pow(t, 2); //Math.pow(x[i], 3) * t + t; // точное решение
                Yn[i] = Math.pow(a, h);//0; /*Math.pow(a, h);*/ //нулевой слой??? но это не точно!

                gamma = a * (t / Math.pow(h, 2));
                A[i] = -gamma;
                C[i] = -(1 + 2 * gamma);
                B[i] = -gamma;
                F[i] = -Yn[i] - t * (2 * x[i] * t - 6 * x[i]); //(Math.pow(x[i], 3) - 6 * x[i] * t + 1);

            }
            mu[1] = 0;
            mu[2] = Math.pow(a, 3) * t + t;
            lam[1] = 0;
            lam[2] = 0;
            alfa[1] = lam[1];
            beta[1] = mu[1];
            for (int i = 1; i <= m - 1; i++) {
                z = C[i] - A[i] * alfa[i];
                alfa[i + 1] = B[i] / z;
                beta[i + 1] = (A[i] * beta[i] + F[i]) / z;
            }
            u[m] = (mu[2] + lam[2] * beta[m]) / (1 - alfa[m] * lam[2]);

            for (int i = m; i >= 1; i--) {
                u[i - 1] = alfa[i] * u[i] + beta[i];
            }

            for (int i = 0; i <= m; i++) {
                System.out.println("Точное: " + u_x[i] + "\t\t" + "Приближенное: " + "  u=" + u[i]);
            }
            /* A[0] = 0;
             B[n - 1] = 0;
             for (int i = 0; i < n; i++) {
             m = A[i] / C[i - 1];
             C[i] = C[i] - m * B[i - 1];
             F[i] = F[i] - m * F[i - 1];
             }
             x[n - 1] = F[n - 1] / C[n - 1];

             for (int i = n - 2; i >= 0; i--) {
             x[i] = (F[i] - B[i] * x[i + 1]) / C[i];
             }*/
        }
    }

    public static void main(String[] args) {
        logic();
    }

}
