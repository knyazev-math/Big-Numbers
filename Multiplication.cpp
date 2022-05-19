#include <complex>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
using namespace std;

int const logLimit = 19;
int const limit = 1 << logLimit;

vector <int> rev;

void calcRev() {
	rev = vector <int>(limit, 0);
	for (int i = 0; i < limit; i++) {
		for (int k = 0; k < logLimit; k++) {
			if (i & (1 << k)) {
				rev[i] ^= 1 << (logLimit - k - 1);
			}
		}
	}
}

using Num = complex <double>;

double const Pi = acos(-1.0);

vector <Num> z;

void calcZ() {
	z = vector <Num>(limit);
	for (int i = 0; i < limit; i++) {
		z[i] = Num(cos(i * 2 * Pi / limit), sin(i * 2 * Pi / limit));
	}
}

vector <Num> fft(const vector <Num>& a0, bool inv) {
	vector <Num> a = a0;
	for (int i = 0; i < limit; i++)
		if (i < rev[i])
			swap(a[i], a[rev[i]]);
	if (inv)
		reverse(z.begin() + 1, z.end());
	for (int k = 0, span = 1, step = limit / 2; k < logLimit;
		k++, span *= 2, step /= 2) {
		for (int i = 0; i < limit; i += 2 * span)
			for (int j = 0; j < span; j++) {
				int u = i + j;
				int v = i + j + span;
				Num x = a[u] + a[v] * z[j * step];
				Num y = a[u] + a[v] * z[j * step + limit / 2];
				a[u] = x;
				a[v] = y;
			}
	}
	if (inv) {
		reverse(z.begin() + 1, z.end());
		for (int i = 0; i < limit; i++)
			a[i] /= limit;
	}
	return a;
}


pair <vector <Num>, pair<int, int>> readNumber() {
	string s;
	cin >> s;
	int zeros = 0, minus = 1, start = 0;
	vector <Num> res(limit, Num(0));
	int n = int(s.size());
	bool zer = true;
	if (s[0] == '-') {
		minus = -1;
		start = 1;
	}
	for (int i = n - 1; i >= start; i--) {
		if (zer) {
			if (s[i] == '0') {
				zeros += 1;
			}
			else {
				zer = false;
			}
		}
	}
	for (int i = start; i < n; i++)
		res[i - start] = Num(s[i] - '0');
	return { res, {zeros, minus} };
}


int main() {
	calcRev();
	calcZ();
	auto a = readNumber();
	auto b = readNumber();
	auto fa = fft(a.first, false);
	auto fb = fft(b.first, false);
	auto fc = vector <Num>(limit);
	for (int i = 0; i < limit; i++)
		fc[i] = fa[i] * fb[i];
	auto c = fft(fc, true);
	vector <int> res(limit);
	long long carry = 0;
	int pos = -1;
	bool started = false;
	for (int i = limit - 1; i > -1; i--) {
		if ((abs(c[i].real()) > 0.3) && (not started)) {
			pos = i;
			started = true;
		}
		if (started) {
			long long cur = (long long)(round(c[i].real())) + carry;
			carry = cur / 10;
			res[i] = cur % 10;
		}
	}
	if ((a.second.second * b.second.second < 0) && (started)) {
		cout << "-";
	}
	if (carry > 0) {
		cout << carry;
	}
	for (int i = 0; i < pos + 1; i++) {
		cout << res[i];
	}
	if (pos == -1) {
		cout << "0";
	}
	else {
		for (int i = 0; i < a.second.first + b.second.first; i++) {
			cout << "0";
		}
	}
	cout << endl;
}
