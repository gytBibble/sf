#include<iostream>
#include<math.h> 
#include<float.h>
#include<stdlib.h>
#include<fstream>
using namespace std;
typedef struct Point {//�������Ͷ��� 
	double x;
	double y;
}point;

typedef struct MinPoints {
	point a;
	point b;
}minpoints;
typedef struct Num_Min_Points {
	double min;
	int num;
}num_min_points;
//����a,b����ľ���
double DstPoints(const point &a, const point &b);
//���η��ݹ���� ��������
void MinAll(point P[], point Q[], int n, minpoints *mp, num_min_points *nmp);
//�鲢����,mark=0ʱ��x��������mark=1ʱ��y��������
void MergeSort(point *A, int n, const int mark);
void Merge(point *A, point *L, int Llength, point *R, int Rlength, const int mark);

int main(void)
{
	int n = 0;
	cout << "����������Ҫ����ĵ�ĸ���\n";
	cin >> n;
	point *test = (point *)malloc(n * sizeof(point));//���Ե㼯 
	point *P = (point *)malloc(n * sizeof(point));//P��
	point *Q = (point *)malloc(n * sizeof(point));//Q��
	minpoints *mp = (minpoints *)malloc(2 * n * sizeof(minpoints));//�洢������
	num_min_points *nmp = (num_min_points *)malloc(sizeof(num_min_points));//��ʾ��ͬ���������Ը�������ʼΪ0
	nmp->min = DBL_MAX;
	nmp->num = 0;
	printf("��������������Ҫ����ĵ㣬�����ʽ(�Ⱥ������������)�磺3.1 5;\n");
	for (int i = 0; i < n; i++)
	{
		cin >> test[i].x >> test[i].y;
		P[i].x = Q[i].x = test[i].x;
		P[i].y = Q[i].y = test[i].y;
	}
	MergeSort(P, n, 0);
	MergeSort(Q, n, 1);
	MinAll(P, Q, n, mp, nmp);
	if (nmp->num)
	{
		cout << "�ҵ�������,����Ϊ��" << nmp->min << endl
			<< "���Ϊ��\n";
		for (int i = 0; i < nmp->num; i++)
		{
			cout << "(" << mp[i].a.x << "," << mp[i].a.y << ") "
				<< "(" << mp[i].b.x << "," << mp[i].b.y << ")\n";
		}
	}
	else {
		cout << "δ�ҵ������ԣ�\n";
	}
	return 0;
}
//����a,b����ľ���
double DstPoints(const point &a, const point &b)
{
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}
//���η��ݹ���� ��������
void MinAll(point P[], point Q[], int n, minpoints *mp, num_min_points *nmp)
{
	double min = DBL_MAX;//��ʼ��Ϊ����� 
	if (n == 1)//ֻ��һ���� 
	{
		return;
	}
	else if (n == 2)//ֻ�������㣬ֱ�ӷ���������룬�޸ĵ������ 
	{
		min = DstPoints(P[0], P[1]);
		if (min < nmp->min)
		{
			mp[0].a = P[0];
			mp[0].b = P[1];
			nmp->num = 1;
			nmp->min = min;
		}
		else if (min == nmp->min)
		{
			mp[nmp->num].a = P[0];
			mp[nmp->num].b = P[1];
			nmp->num += 1;
		}
		return;
	}
	else
	{
		int m = (n + 1) / 2;
		point *PL = (point *)malloc(m * sizeof(point));
		point *PR = (point *)malloc((n - m) * sizeof(point));
		point *QL = (point *)malloc(n * sizeof(point));
		point *QR = (point *)malloc(n * sizeof(point));
		//��P��Ϊ�����֣�PL��PR
		for (int i = 0; i < m; i++)
		{
			PL[i] = P[i];
		}
		for (int i = 0; i < n - m; i++)
		{
			PR[i] = P[m + i];
		}
		//����Q��QL��QR��
		for (int i = 0; i < n; i++)
		{
			QL[i] = QR[i] = Q[i];
		}
		//ɾ��QL��QR�и��Բ��ڷ�Χ�ڵĵ�
		for (int i = 0, j = 0; i < n; i++)
		{
			if (QL[i].x >= PL[0].x && QL[i].x <= PL[m - 1].x)
				QL[j++] = QL[i];
		}
		for (int i = 0, j = 0; i < n; i++)
		{
			if ((QR[i].x > PR[0].x && QR[i].x <= PR[n - m - 1].x) || (QR[i].x == PR[0].x&&QR[i].y >= PR[0].y))
				QR[j++] = QR[i];
		}

		//��PL��PR��QL��QR����ݹ���̴���
		MinAll(PL, QL, m, mp, nmp);
		MinAll(PR, QR, n - m, mp, nmp);
		
		//ɨ�豾��Q��,ɾ����x���겻�ڴ��ڵ����е�
		int strip_num = 0;
		for (int i = 0; i < n; i++)
		{
			if (fabs(Q[i].x - P[m - 1].x) <= min)
				Q[strip_num++] = Q[i];
		}
		double minflag = min;
		double minLR = DBL_MAX;
		//��Q�����µĵ���������
		for (int i = 0; i < strip_num; i++)
		{
			for (int j = i + 1; j < strip_num && (Q[j].y - Q[i].y) <= minflag; j++)//������7�Σ���ΪO(1)
			{
				minLR = DstPoints(Q[i], Q[j]);
				if (minLR < min)//С�ڣ��޸���ĿΪ1
				{
					min = minLR;
					mp[0].a = Q[i];
					mp[0].b = Q[j];
					nmp->min = min;
					nmp->num = 1;
				}
				else if (minLR == min)
				{
					mp[nmp->num].a = Q[i];
					mp[nmp->num].b = Q[j];
					nmp->num += 1;
				}
			}
		}
		return;
	}
}
//�鲢����,mark=0ʱ��x��������mark=1ʱ��y��������
void MergeSort(point *A, int n, const int mark)
{
	if (n < 2) return;
	int mid = n / 2;
	point *L = new point[mid];
	point *R = new point[n - mid];
	for (int i = 0; i < mid; i++) L[i] = A[i];
	for (int i = mid; i < n; i++) R[i - mid] = A[i];
	MergeSort(L, mid, mark);
	MergeSort(R, n - mid, mark);
	Merge(A, L, mid, R, n - mid, mark);
	delete[]L;
	delete[]R;
}
void Merge(point *A, point *L, int Llength, point *R, int Rlength, const int mark)
{
	int i, j, k;
	i = 0; j = 0; k = 0;
	while (i<Llength && j< Rlength) {
		if (mark == 0 && (L[i].x < R[j].x || (L[i].x == R[j].x && L[i].y < R[j].y))) A[k++] = L[i++];
		else if (mark == 1 && L[i].y < R[j].y) A[k++] = L[i++];
		else A[k++] = R[j++];
	}
	while (i < Llength) A[k++] = L[i++];
	while (j < Rlength) A[k++] = R[j++];
}
