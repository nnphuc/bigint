#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <assert.h>

/**************************
 *
 *  can use:
 * add(a,b,a)
 * sub(a,b,a)
 * div_ll(a,x,a)
 * mul_fft(a,b,a)
 * _mul(a,b,a)
 * mul(a,b,a)
 * divmod(a,b,a,r)
 * recrusive_divmod(a,b,q,r)
 *
 *****************************/
typedef long long ll;
#define  real long double
typedef real complex cplx;
real PI;
cplx *memo;
int FFTN;
void init_fft(int n){
	PI=(real)atan2l(1,1)*4;
	int i;
	memo=malloc(sizeof(cplx)*n);
	FFTN=n;
	for(i=0;i<n;++i){
		memo[i]=cexpl(-I*PI*i/FFTN);
	}
}
void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 		int i;
		for ( i = 0; i < n; i += 2 * step) {
			cplx t =memo[i*FFTN/n] * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
void _ifft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) {
		_ifft(out, buf, n, step * 2);
		_ifft(out + step, buf + step, n, step * 2);
 		int i;
		for ( i = 0; i < n; i += 2 * step) {
			cplx t =conjl(memo[i*FFTN/n]) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
void fft(cplx buf[], int n)
{
	cplx out[n];
	int i;
	for ( i = 0; i < n; i++) out[i] = buf[i];
	_fft(buf, out, n, 1);
}

void ifft(cplx buf[], int n)
{
	cplx out[n];
	int i;
	for ( i = 0; i < n; i++) out[i] = buf[i];
	_ifft(buf, out, n, 1);
	for ( i = 0; i < n; i++) buf[i] =buf[i]/n;
}

void show_cplx(char *s,cplx *a,int n){
	int i;
	puts(s);
	for(i=0;i<n;++i){
		if (!cimagl(a[i]))
			printf("%Lf ", creall(a[i]));
		else
			printf("(%Lf, %Lf) ", creall(a[i]), cimagl(a[i]));
	}
	puts("");
}


#define B 1000000
#define maxn 2000
#define C 6
#define FIX "6"
#define max2(a,b) ((a)>(b)?a:b)
#define min2(a,b) ((a)<(b)?a:b)
#define CMP(a,b) ((a)<(b)?-1:(a)>(b)?1:0)

typedef struct num{
	ll a[maxn];
	int sign;
	int n;
} num;

void trim(num *n){
	while (n->n && n->a[n->n-1]==0){
		--n->n;
	}
	if(!n->n){
		n->sign=1;
	}
}
void print(const num* n){
	int i=n->n-1;
	if(n->n==0){
		printf("0\n");
		return;
	}
	if(n->sign<0) printf("-");
	printf("%lld",n->a[i]);
	while(--i>=0) printf("%0"FIX"lld",n->a[i]);
	printf("\n");
}
ll s_to_ll(char *s,int st,int nd){
	int i;
	ll r=0;
	for(i=st;i<=nd;++i) r=r*10+(ll)(s[i]-'0');
	return r;
}
void s_to_num(char *s,num *n){
	char *p=s;
	if(*p=='-'){
		n->sign=-1;
		p++;
	}else n->sign=1;
	int len=strlen(p);
	int nd=len-1;
	n->n=0;
	while(nd>=0){
		n->a[n->n++]=s_to_ll(p,(nd-C+1>=0?nd-C+1:0),nd);
		nd-=C;
	}
}

void ll_to_num(ll x, num *n){
	if(x<0){
		n->sign=-1;
		x=-x;
	}else n->sign=1;
	n->n=0;
	while (x){
		n->a[n->n++]=x%B;
		x/=B;
	}
}
void rand_num(num *n,int len){
	//random num with len digit;
	srand(time(NULL));
	char *s=malloc(len+2);
	int i;
	s[0]=rand()%9+1+'0';
	for(i=1;i<len;++i){
		s[i]=rand()%10+'0';
	}
	s[len]='\0';
	s_to_num(s,n);
	free(s);
}

void _add(num *a,num *b,num *c){
	int i;

	ll m=0;
	for(i=0;i<max2(a->n,b->n);++i){
		m+= (i<a->n?a->a[i]:0) + (i<b->n?b->a[i]:0) ;
		c->a[i]=m%B;
		m/=B;
	}
	c->n=max2(a->n,b->n);
	while(m){
		c->a[c->n++]=m%B;
		m/=B;
	}
}

void shift_right(num *a, int n){
	int i;
	if(!a->n) return;
	for(i=a->n-1;i>=0;--i){
		a->a[i+n]=a->a[i];
	}
	for(i=0;i<n;++i){
		a->a[i]=0;
	}
	a->n+=n;
}

void _sub(num *a,num *b ,num *c){
	int i;
	ll t=0,m=0;
	for(i=0;i<a->n;++i){
		t=a->a[i] - (i<b->n?b->a[i]:0) +m;
		if(t<0){
			c->a[i]=t+B;
			m=-1;
		}else{
			c->a[i]=t;
			m=0;
		}
	}
	c->n=a->n;
	trim(c);
}

void mul_ll(num *a,ll x,num *c){
	if(x==0){
		c->n=0;
		c->sign=1;
		return;
	}
	c->sign=a->sign;
	if(x<0){
		c->sign=-c->sign;
		x=-x;
	}
	int i;
	ll t=0;
	for(i=0;i< a->n;++i){
		t += a->a[i]*x;
		c->a[i]=t%B;
		t/=B;
	}
	c->n=a->n;
	while(t){
		c->a[c->n++]=t%B;
		t/=B;
	}
	trim(c);
}

void copy(num *dst,num *src){
	dst->n=src->n;
	dst->sign=src->sign;
	int i;
	for(i=0;i<dst->n;++i) dst->a[i]=src->a[i];
}

void _mul(num *a, num *b,num *c1){
	num c;
	c.n=0;
	num t,tt;
	int i;
	for(i=b->n-1;i>=0;--i){
		mul_ll(a,b->a[i],&t);
		shift_right(&c,1);
		_add(&c,&t,&c);
	}
	c.sign=a->sign*b->sign;
	copy(c1,&c);
}
void mul_fft(num *a,num *b,num *c){
	int i,n1=a->n,n2=b->n,n=1;
	while(n<a->n+b->n){
		n<<=1;
	}
	cplx fa[n],fb[n];
	for(i=0;i<n1;++i){
		fa[i]=a->a[i];
	}
	for(;i<n;++i){
		fa[i]=0;
	}
	for(i=0;i<n2;++i){
		fb[i]=b->a[i];

	}
	for(;i<n;++i){
		fb[i]=0;
	}

	fft(fa,n);
	fft(fb,n);

	for(i=0;i<n;++i)fa[i]*=fb[i];

	ifft(fa,n);
	ll m=0;
	//show_cplx("res",fa,n);
	for(i=0;i<n;++i){
		m+=(ll)roundl(creall(fa[i]));
		c->a[i]=m%B;
		m/=B;
	}
	c->n=n;
	while(m){
		c->a[c->n++]=m%B;
		m/=B;
	}
	trim(c);
	c->sign=a->sign*b->sign;
}
void mul(num* a,num *b,num *c){
	if(a->n<3 ||b->n<3){
		_mul(a,b,c);
	}else{
		mul_fft(a,b,c);
	}
}
int _cmp(num *a,num *b){
	if(a->n!=b->n) return CMP(a->n,b->n);
	int i;
	for(i=a->n-1;i>=0;--i){
		if(a->a[i]!=b->a[i]) return CMP(a->a[i],b->a[i]);
	}
	return 0;
}
void add(num *a,num *b,num *c){
	if(a->sign==b->sign){
		_add(a,b,c);
		c->sign=a->sign;
	}else{
		int t=_cmp(a,b);
		if(t==0){
			ll_to_num(0,c);
		}else if(t>0){
			_sub(a,b,c);
			c->sign=a->sign;
		}else{
			_sub(b,a,c);
			c->sign=b->sign;
		}
	}
}
void sub(num *a,num *b,num *c){
	if(a->sign!=b->sign){
		_add(a,b,c);
		c->sign=a->sign;
	}else{
		int t=_cmp(a,b);
		if(t==0){
			ll_to_num(0,c);
		}else if(t>0){
			_sub(a,b,c);
			c->sign=a->sign;
		}else{
			_sub(b,a,c);
			c->sign=-a->sign;
		}
	}
}

num ONE;
ll div_ll(num *a,ll x,num *c){
	c->sign=a->sign;
	if(x<0){
		c->sign=-c->sign;
		x=-x;
	}
	int i;
	ll cur,rem=0;
	for(i=a->n-1;i>=0;--i){
		cur=a->a[i] + rem*B;
		c->a[i]=cur/x;
		rem=cur%x;
	}
	c->n=a->n;
	trim(c);
	return rem;
}
void divmod(num *a1,num *b1,num *q,num *r){
	if(b1->n<=2){
		ll tmp=b1->sign*((b1->n>0?b1->a[0]:0)+(b1->n==2?b1->a[1]:0)*B);
		ll rem=div_ll(a1,tmp,q);
		ll_to_num(rem,r);
		return;
	}
	ll norm=B/(b1->a[b1->n-1]+1);
	num a,b;
	mul_ll(a1,norm,&a);
	mul_ll(b1,norm,&b);
	a.sign=1;
	b.sign=1;
	int i;
	num rc;
	r->n=0;
	num t1;
	q->n=a.n;

	for(i=a.n-1;i>=0;--i){
		mul_ll(r,B,&rc);
		ll_to_num(a.a[i],&t1);
		add(&rc,&t1,r);
		ll s1 = r->n <= b.n ? 0 : r->a[b.n];
		ll s2 = r->n <= b.n - 1 ? 0 : r->a[b.n - 1];
		ll d=(B*s1+s2)/b.a[b.n-1];
		mul_ll(&b,d,&t1);
		sub(r,&t1,r);
		while(!r->n ||(r->sign<0 )||r->a[r->n-1]==0){
			add(r,&b,r);
			d--;
		}
		q->a[i]=d;
	}
	trim(q);
	div_ll(r,norm,r);
	q->sign=a1->sign*b1->sign;
	r->sign=a1->sign;
	trim(r);

}
void recursive_divmod(num *a,num *b,num *q,num *r){
	if(b->n<=2){
		ll tmp=b->sign*((b->n>0?b->a[0]:0)+(b->n==2?b->a[1]:0)*B);
		div_ll(a,tmp,q);
		num tmp_mul;
		mul_ll(q,tmp,&tmp_mul);
		sub(a,&tmp_mul,r);
		return;
	}
	int i,k=b->n/2;
	num a1,a0,b1,b0,q1,r1;
	num ta,t;
	b1.n=0;
	b1.sign=b->sign;

	for(i=k;i<b->n;++i){
		b1.a[b1.n++]=b->a[i];
	}
	b0.n=0;
	b0.sign=b->sign;
	for(i=0;i<k;++i){
		b0.a[b0.n++]=b->a[i];
	}
	ta.n=0;
	ta.sign=a->sign;
	for(i=2*k;i<a->n;++i){
		ta.a[ta.n++]=a->a[i];
	}
	a0.n=0;
	a0.sign=a->sign;
	for(i=0;i<min2(2*k,a->n);++i){
		a0.a[a0.n++]=a->a[i];
	}
	recursive_divmod(&ta,&b1,&q1,&r1);

	shift_right(&r1,2*k);
	add(&r1,&a0,&t);
	copy(&ta,&t);
	mul(&q1,&b0,&t);
	shift_right(&t,k);
	num tb;
	sub(&ta,&t,&tb);
	copy(&ta,&tb);
	copy(&tb,b);
	shift_right(&tb,k);


	while(ta.sign<0){
		sub(&q1,&ONE,&t);
		copy(&q1,&t);
		add(&ta,&tb,&t);
		copy(&ta,&t);
	}

	num tta;
	tta.n=0;
	tta.sign=ta.sign;
	for(i=k;i<ta.n;++i){
		tta.a[tta.n++]=ta.a[i];
	}

	num q0,r0;
	a0.n=0;
	a0.sign=ta.sign;
	for(i=0;i<min2(k,ta.n);++i){
		a0.a[a0.n++]=ta.a[i];
	}
	recursive_divmod(&tta,&b1,&q0,&r0);

	shift_right(&r0,k);
	add(&a0,&r0,&t);

	copy(&tta,&t);
	mul(&q0,&b0,&t);
	sub(&tta,&t,&a0);
	copy(&tta,&a0);

	while(tta.sign<0){
		sub(&q0,&ONE,&t);
		copy(&q0,&t);
		add(&tta,b,&t);
		copy(&tta,&t);

	}
	shift_right(&q1,k);
	add(&q1,&q0,q);
	copy(r,&tta);
}
void sqrt_num(num *a,num *c){

	num x;
	//ll_to_num(1,&x);
	rand_num(&x,a->n/2*C+1);
	num y,yc,t,tt;
	for(;;){
		recursive_divmod(a,&x,&yc,&tt); // yc=a/x
		//puts("---");
		//print(&yc);
		add(&yc,&x,&y); //y=a/x+x
		//print(&x);
		//print(&y);

		div_ll(&y,2,&yc); // yc=(a/x+x)/2
		sub(&x,&yc,&t);
		//print(&yc);
		if(_cmp(&t,&ONE)==0 && t.sign<0 ||t.n==0){
			copy(c,&x);
			return;
		}
		copy(&x,&yc);
	}
}


int main(void) {
	init_fft(1<<12);
	ll_to_num(1,&ONE);
	char s[1000];
	scanf("%s",s);
	int i,j;
	num a,b,c,x;
	s_to_num(s,&a);
	scanf("%s",s);
	s_to_num(s,&b);
	recursive_divmod(&a,&b,&a,&x);
	print(&a);
	print(&x);
	return 0;
}


