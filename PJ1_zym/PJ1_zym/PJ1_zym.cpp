#include <iostream>
#include <stdio.h>
#include <ctime>
#include <cmath>
#define dens 1000
#define eps 1e-8
using namespace std;
void inputPolygon();
void precomputation();
void inputRecords();
void func(int oid,double x,double y,int t);
int func1(double x,double y);
void make_grid(double xl,double xr,double yd,double yu,int a,int b);
void outputResults();
double px[202],py[202],rx[20000000],ry[20000000],gx[dens+1],gy[dens+1];
int oid,t,M,N,grid[dens][dens];
double maxx=-10000,maxy=-10000,minx=10000,miny=10000,vmax,kiv[102],k[102];
bool ans[20000000];
int func1(double x,double y)
{
	int count = 0;
	if (x<=minx || x>=maxx || y<=miny || y>=maxy)
		return 0;  //外接矩形之外
	for (int i=0;i<N;i++)
	{
		if (y < py[i] && y < py[(i+1)]) continue;
		if (y > py[i] && y > py[(i+1)]) continue;
		if (x > px[i] && x > px[(i+1)]) continue;
		if (y==py[i])
		{
			if (x==px[i])
				return 0;
			if (y!=py[(i+1)])
			{
				if ((py[(i+N-1)]>y && py[(i+1)]>y) || (py[(i+N-1)]<y && py[(i+1)]<y))
					count++;
				continue;
			}
			else
			{
				if ((x > px[i] && x < px[(i+1)]) || (x < px[i] && x > px[(i+1)]))
					return 0;
				i++;
				if ((py[(i+N-2)]>y && py[(i+1)]>y) || (py[(i+N-2)]<y && py[(i+1)]<y))
					count++;
				continue;
			}
		}
		double ix = (y-py[i]) * kiv[i] +px[i];
		//cout << x << " " << ix << endl;
		if (fabs(ix-x)<eps)
			return 0;
		if (ix>x) count++;
	}
	return count%2;
}
void make_grid(double xl,double xr,double yd,double yu,int a,int b)
{
	for (int i=0;i<N;i++)
	{
		if (yd > py[i] && yd > py[(i+1)]) continue;
		if (yu < py[i] && yu < py[(i+1)]) continue;
		if (xl > px[i] && xl > px[(i+1)]) continue;
		if (xr < px[i] && xr < px[(i+1)]) continue;
		double ix1 = (yd-py[i]) * kiv[i] +px[i];
		double ix2 = (yu-py[i]) * kiv[i] +px[i];
		double iy1 = (xl-px[i]) * k[i] +py[i];
		double iy2 = (xr-px[i]) * k[i] +py[i];
		if (ix1-xl>eps && xr-ix1>eps)
		{
			grid[a][b]=2;
			return;
		}
		if (ix2-xl>eps && xr-ix2>eps)
		{
			grid[a][b]=2;
			return;
		}
		if (iy1-yd>eps && yu-iy1>eps)
		{
			grid[a][b]=2;
			return;
		}
		if (iy1-yd>eps && yu-iy1>eps)
		{
			grid[a][b]=2;
			return;
		}
	}
	grid[a][b]=func1((xl+xr)/2,(yd+yu)/2);
	return;
}
void inputPolygon()
{
	freopen("Polygon.txt","r",stdin);
	scanf("%d",&N);
	for (int i=0;i<N;i++)
		scanf("%lf%lf",&px[i],&py[i]);
	fclose(stdin);
}
void precomputation()
{ 
	for (int i=0;i<N;i++)
	{
		px[i+N]=px[i];
		py[i+N]=py[i];
		if (px[i] > maxx) maxx=px[i];
		if (px[i] < minx) minx=px[i];
		if (py[i] > maxy) maxy=py[i];
		if (py[i] < miny) miny=py[i];
		if (py[(i+1)]!=py[i])
			kiv[i]=(px[(i+1)]-px[i]) /(py[(i+1)]-py[i]);
		else kiv[i]=1e10;
		if (px[(i+1)]!=px[i])
			k[i]=(py[(i+1)]-py[i]) /(px[(i+1)]-px[i]);
		else k[i]=1e10;
	}
	for (int i=0;i<dens;i++)
	{
		gx[i]=minx+(maxx-minx)*i/dens;
		gy[i]=miny+(maxy-miny)*i/dens;
	}
	gx[dens]=maxx;
	gy[dens]=maxy;
	for (int i=0;i<dens;i++)
		for (int j=0;j<dens;j++)
			make_grid(gx[i],gx[i+1],gy[j],gy[j+1],i,j);
}
void inputRecords()
{
	freopen("Records.txt","r",stdin);
	scanf("%lf%d",&vmax,&M);
	for (int i=0;i<M;i++)
		scanf("%d%lf%lf%d",&oid,&rx[i],&ry[i],&t);
	fclose(stdin);
}
void func(int oid,double x,double y,int t)
{
	if (x<=minx || x>=maxx || y<=miny || y>=maxy)
	{
		ans[oid]=0;
		return;  //外接矩形之外
	}
	if (grid[(int) ((x-minx)*dens/(maxx-minx))][(int) ((y-miny)*dens/(maxy-miny))]<2)
	{
		ans[oid]=(grid[(int) ((x-minx)*dens/(maxx-minx))][(int) ((y-miny)*dens/(maxy-miny))]==1);
		return;
	}
	int count = 0;
	for (int i=0;i<N;i++)
	{
		if (y < py[i] && y < py[(i+1)]) continue;
		if (y > py[i] && y > py[(i+1)]) continue;
		if (x > px[i] && x > px[(i+1)]) continue;
		if (y==py[i])
		{
			if (x==px[i])
			{
				ans[oid]=0;
				return;
			}
			if (y!=py[(i+1)])
			{
				if ((py[(i+N-1)]>y && py[(i+1)]>y) || (py[(i+N-1)]<y && py[(i+1)]<y))
					count++;
				continue;
			}
			else
			{
				if ((x > px[i] && x < px[(i+1)]) || (x < px[i] && x > px[(i+1)]))
				{
					ans[oid]=0;
					return;
				}
				i++;
				if ((py[(i+N-2)]>y && py[(i+1)]>y) || (py[(i+N-2)]<y && py[(i+1)]<y))
					count++;
				continue;
			}
		}
		double ix = (y-py[i]) * kiv[i] +px[i];
		//cout << x << " " << ix << endl;
		if (fabs(ix-x)<eps)
		{
			ans[i]=0;
			return;
		}
		if (ix>x) count++;
	}
	ans[oid]=(count%2==1);
}
void outputResults()
{
	freopen("Result.txt","w",stdout);
	for (int i=0;i<M;i++)
	{
		if (ans[i]) printf("true\n");
		else printf("false\n");
	}
	fclose(stdout);
}
int main()
{
	inputPolygon(); //读入多边形
	precomputation(); //对二维平面和多边形建立索引
	inputRecords(); //读入点记录 
	for (int i=0; i<M; ++i)
	{
		func(i,rx[i],ry[i], t); //核心函数，判断点记录对应的点当时是否位于多边形内
	}
	outputResults(); //输出结果	
	return 0;
}