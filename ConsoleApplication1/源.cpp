#include <iostream>

#define ROWS 512
#define COLS 512
#define PI 3.14159265

using namespace std;

int occurs_count[64][10000];
double adct[ROWS][COLS];


void dct_trans()
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
					inblock[u][v] = inblock[u][v] - 128;//£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡

			/* Forward DCT of one block */
			double temp = 0.0;
			for (int u = 0; u < 8; u++)
				for (int v = 0; v < 8; v++)
				{
					temp = 0.0;
					for (int x = 0; x < 8; x++)
						for (int y = 0; y < 8; y++)
							temp += (double)inblock[x][y] * cos((2 * x + 1)*u*PI/ 16)*cos((2 * y + 1)*v*PI / 16);
					adct[i * 8 + u][j * 8 + v] = temp*c[u] * c[v] / 4;
				}
		}
	}
}

int main()
{
	dct_trans();
	return 0;
}