#pragma once
#include <iostream>
#include <cfloat>

#define ROWS 512
#define COLS 512
#define PI 3.14159265
#define upper(x) (((x) + 0.5) / 1)
#define lower(x) ((x) / 1)
#define sqr(x) ((x) * (x))
#define QMAX 150
#define BPPSCALE 1000
#define INFINITY DBL_MAX

using namespace std;

short int zigzag[8][8] = { 0, 1, 5, 6,14,15,27,28,
2, 4, 7,13,16,26,29,42,
3, 8,12,17,25,30,41,43,
9,11,18,24,31,40,44,53,
10,19,23,32,39,45,52,54,
20,22,33,38,46,51,55,60,
21,34,37,47,50,56,59,61,
35,36,48,49,57,58,62,63 };
static int occurs_count[64][10000] = { { 0 } };
double adct[ROWS][COLS];
int R[64][QMAX];
double D[64][QMAX];
int v_max, v_min;
int max_rate = 0;
double *least_D[64];
int *choice[64];

void dctTrans()
{
	FILE *infd;
	unsigned char input[ROWS][COLS];
	int inblock[8][8];
	double c[8] = { 0.707,1.0,1.0,1.0,1.0,1.0,1.0,1.0 };
	/* input the original pgm image to input[][] */
	infd = fopen("lena512.raw", "rb");
	fread(input, sizeof(unsigned char), ROWS*COLS, infd);
	fclose(infd);

	/* the ROWSxCOLS image is processed by ROWS/8xCOLS/8 block lines and columns */
	for (int i = 0; i < ROWS / 8; i++)
	{
		for (int j = 0; j < COLS / 8; j++)
		{/*for each image block, go through following step */

		 /* change the format of unsigned char to integer before computation*/
			for (int u = 0; u < 8; u++)
				for (int v = 0; v < 8; v++)
					inblock[u][v] = (int)input[i * 8 + u][j * 8 + v];

			/* remove the mean of pixels in the image block*/
			for (int u = 0; u < 8; u++)
				for (int v = 0; v < 8; v++)
					inblock[u][v] = inblock[u][v] - 128;

			/* Forward DCT of one block */
			double temp = 0.0;
			for (int u = 0; u < 8; u++)
				for (int v = 0; v < 8; v++)
				{
					temp = 0.0;
					for (int x = 0; x < 8; x++)
						for (int y = 0; y < 8; y++)
							temp += (double)inblock[x][y] * cos((2 * x + 1)*u*PI / 16)*cos((2 * y + 1)*v*PI / 16);
					adct[i * 8 + u][j * 8 + v] = temp*c[u] * c[v] / 4;
					if (v_max < adct[i * 8 + u][j * 8 + v])
						v_max = adct[i * 8 + u][j * 8 + v];
					if (v_min > adct[i * 8 + u][j * 8 + v])
						v_min = adct[i * 8 + u][j * 8 + v];
				}
		}
	}
}

void getherStats()
{
	for (int i = 0; i < ROWS / 8; i++)
	{
		for (int j = 0; j < COLS / 8; j++)
		{
			int v;
			for (int x = 0; x < 8; x++)
			{
				for (int y = 0; y < 8; y++)
				{
					if (adct[i * 8 + x][j * 8 + y] >= 0)
						v = lower(2 * adct[i * 8 + x][j * 8 + y]);
					else
						v = -1 * lower(-2 * adct[i * 8 + x][j * 8 + y]);
					occurs_count[zigzag[x][y]][v + 5000]++;//加5000，把0移到5000的位置，小于5000表示负数
				}
			}
		}
	}
}

void fillR()
{
	int F = (int)ROWS*COLS / 64;
	for (int x = 0; x < 8; x++)
	{
		for (int y = 0; y < 8; y++)
		{
			for (int q = 1; q < QMAX; q++)
			{
				double entropy = 0;
				for (int i = -(int)(-v_min / (double)q + 0.5); i <= (int)(v_max / (double)q + 0.5); i++)//遍历每一个可能的量化后的系数
				{
					int count = 0;
					for (int v = (int)(v_min * 2) - 2; v < (int)(v_max * 2) + 2; v++)//遍历每一个v，范围是-2VMAX到2VMAX
					{
						if (v < 0 && -(int)(-v / (2.*q) + 0.5) == i)
						{
							count += occurs_count[zigzag[x][y]][v + 5000];
						}
						else if (v >= 0 && (int)(v / (2.*q) + 0.5) == i)
						{
							count += occurs_count[zigzag[x][y]][v + 5000];
						}
					}
					double prob = count / (double)F;
					if (prob > 0)
					{
						entropy = entropy - (prob * (log(prob) / log(2)));
					}
				}

				R[zigzag[x][y]][q] = (int)(entropy / 64. * BPPSCALE);
				if (R[zigzag[x][y]][q] > max_rate)
					max_rate = R[zigzag[x][y]][q];
			}
		}
	}
}

void fillD()
{
	for (int x = 0; x < 8; x++)
	{
		for (int y = 0; y < 8; y++)
		{
			for (int q = 1; q < QMAX; q++)
			{
				D[zigzag[x][y]][q] = 0;
				double original_val;
				int quantized_val;
				double error;
				for (int v = (int)(v_min * 2) - 2; v < (int)(v_max * 2) + 2; v++)//遍历每一个v，范围是-2VMAX到2VMAX
				{
					original_val = v / 2. + ((v < 0) ? -0.25 : 0.25);
					quantized_val = (v < 0) ? -(int)(-v / (2.*q) + 0.5) : (int)(v / (2.*q) + 0.5);
					error = occurs_count[zigzag[x][y]][v + 5000] * sqr(original_val - q * quantized_val);
					D[zigzag[x][y]][q] += error;
				}
				//D[zigzag[x][y]][q] /= N
			}
		}
	}
}