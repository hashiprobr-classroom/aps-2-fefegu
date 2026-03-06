#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}

void fft(complex s[], complex t[], int n, int sign) {

    if (n == 1) {
        t[0] = s[0];
        return;
    }

    complex sp[n/2];
    complex si[n/2];

    for (int i = 0; i < n / 2; i++) {
        sp[i] = s[2 * i]; 
        si[i] = s[2 * i + 1];
    }

    complex tp[n/2];
    complex ti[n/2];

    fft(sp,tp,n/2,sign);
    fft(si,ti,n/2,sign);

    for (int k = 0; k < n / 2; k++) {
        double x = sign * 2 * PI * k / n;

        double cosx = cos(x);
        double sinx = sin(x);

        complex t1;
        t1.a = tp[k].a + cosx * ti[k].a - sinx * ti[k].b;
        t1.b = tp[k].b + cosx * ti[k].b + sinx * ti[k].a;

        complex t2;
        t2.a = tp[k].a - cosx * ti[k].a + sinx * ti[k].b;
        t2.b = tp[k].b - cosx * ti[k].b - sinx * ti[k].a;

        t[k] = t1;
        t[k + n / 2] = t2;
    }
}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    complex c[height];
    complex tc[height];
    for(int i = 0; i < width; i++){ 
        for(int k = 0; k < height;k++){
            c[k] = matrix[k][i];
        }
        fft_forward(c,tc,height);

        for(int k = 0; k < height;k++){
            matrix[k][i] = tc[k];
        }
    }

    complex tl[width];
    for (int i = 0; i < height; i++) { 
        fft_forward(matrix[i], tl, width);

        for (int k = 0; k < width; k++) {
            matrix[i][k] = tl[k];
        }
    }
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    complex tl[width];
    for (int i = 0; i < height; i++) { 
        fft_inverse(matrix[i], tl, width);

        for (int k = 0; k < width; k++) {
            matrix[i][k] = tl[k];
        }
    }

    complex c[height];
    complex tc[height];
    for(int i = 0; i < width; i++){ 
        for(int k = 0; k < height;k++){
            c[k] = matrix[k][i];
        }
        fft_inverse(c,tc,height);

        for(int k = 0; k < height;k++){
            matrix[k][i] = tc[k];
        }
    }
}

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
