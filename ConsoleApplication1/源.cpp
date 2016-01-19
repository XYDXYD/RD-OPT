#include "fun.h"

int main()
{
	int target_rate;

	dctTrans();
	getherStats();
	fillR();
	fillD();

	for (int i = 0; i < 64; i++)
	{
		least_D[i] = new double[max_rate * 64];
		choice[i] = new int[max_rate * 64];
		for (int j = 0; j < max_rate * 64; j++)
		{
			choice[i][j] = 150;
		}
	}
	for (int n = 0; n < 64; n++)//³õÊ¼»¯least_D
	{
		for (int s = 0; s <= max_rate; s++)
		{
			least_D[n][s] = INFINITY;
		}
	}

	for (int q = 1; q < QMAX; q++)
	{
		if (D[0][q] < least_D[0][R[0][q]])
		{
			least_D[0][R[0][q]] = D[0][q];
			choice[0][R[0][q]] = q;
		}
	}
	for (int n = 1; n < 64; n++)
	{
		for (int q = 1; q < QMAX; q++)
		{
			for (int s = 0; s <= max_rate; s++)
			{
				if (D[n][q] + least_D[n - 1][s] < least_D[n][s + R[n][q]])
				{
					least_D[n][s + R[n][q]] = D[n][q] + least_D[n - 1][s];
					choice[n][s + R[n][q]] = q;
				}
			}
		}
	}

	for (int i = 0; i < 64; i++)
	{
		for (int j = 1; j < 64 * max_rate; j++)
		{
			if (choice[i][j] > choice[i][j - 1])
				choice[i][j] = choice[i][j - 1];
			//  cout << choice[i][j] << " ";
		}
		//  cout << endl;
	}

	target_rate = 110;
	for (int n = 63; n >= 0; n--)
	{
		cout << choice[n][target_rate] << "  ";
		target_rate -= R[n][choice[n][target_rate]];
	}

	for (int i = 0; i < 64; i++)
	{
		delete[] least_D[i];
		delete[] choice[i];
	}
	return 0;
}