#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>

using namespace std;

double initial_price = 91.90;
double r = 0.00819347502;
double t = 457.0 / 365;
double q = 0.0479;
double sigma = 0.45116;
double coupon = 0.0325;

int t_step = 100000;
int p_step = 1000;


int main() {
	double dt = t / t_step;
	double ds = 3 * initial_price / p_step;

	vector<int> review_dates;
	review_dates.push_back(92.0 / 457.0 * t_step);
	review_dates.push_back(184.0 / 457.0 * t_step);
	review_dates.push_back(275.0 / 457.0 * t_step);
	review_dates.push_back(365.0 / 457.0 * t_step);

    //grid
    int j0 = p_step / 3;
    int jb = 0.7 * j0;
	vector<vector<double>> vt(t_step, vector<double>(p_step, 0));
	vector<vector<double>> v(t_step, vector<double>(p_step, 0));


	//boundary for vt
	for (int i = 0; i < review_dates[0]; i++) {
		vt[i][0] = 0;
		vt[i][p_step - 1] = (1000 + 1000 * coupon * 1) * exp(-r * (review_dates[0] - i) * dt);
	}
	for (int i = review_dates[0]; i < review_dates[1]; i++) {
		vt[i][0] = 0;
		vt[i][p_step - 1] = (1000 + 1000 * coupon * 2) * exp(-r * (review_dates[1] - i) * dt);
	}
	for (int i = review_dates[1]; i < review_dates[2]; i++) {
		vt[i][0] = 0;
		vt[i][p_step - 1] = (1000 + 1000 * coupon * 3) * exp(-r * (review_dates[2] - i) * dt);
	}
	for (int i = review_dates[2]; i < review_dates[3]; i++) {
		vt[i][0] = 0;
		vt[i][p_step - 1] = (1000 + 1000 * coupon * 4) * exp(-r * (review_dates[3] - i) * dt);
	}
	for (int i = review_dates[3]; i < t_step; i++) {
		vt[i][0] = 0;
		vt[i][p_step - 1] = (1000 + 1000 * coupon * 5) * exp(-r * (t_step - i) * dt);
	}

	//termianl condition for vt
	for (int j = j0; j < p_step; j++) {
		vt[t_step - 1][j] = (1000 + 1000 * coupon * 5) * exp(-r * 5.0 / 365);
	}
	for (int j = jb; j < j0; j++) {
		vt[t_step - 1][j] = (1000 + 1000 * coupon * 5) * exp(-r * 5.0 / 365);
	}
	for (int j = 0; j < jb; j++) {
		vt[t_step - 1][j] = 1000 * (1.0 * j / j0) * exp(-r * 5.0 / 365);
	}

	//terminal condition for v
	for (int j = jb; j < p_step -1; j++) {
		v[t_step - 1][j] = (1000 + 1000 * coupon * 5) * exp(-r * 5.0 / 365);
	}
	for (int j = 0; j < jb; j++) {
		v[t_step - 1][j] = 1000 * (1.0 * j / j0) * exp(-r * 5.0 / 365);
	}


	//solve vt and v backwards
	for (int i = t_step - 2; i >= 0; i--) {
		int coupon_no = 1;
		
        // first review date
		if (i == review_dates[0]) {
			for (int j = 1; j < p_step - 1; j++) {
				vt[i][j] = 0.5 * dt * (sigma * sigma * j - (r - q)) * j * vt[i + 1][j - 1] + 1 - dt * (sigma * sigma * j * j + r) * vt[i + 1][j] + 0.5 * dt * (sigma * sigma * j + (r - q)) * j * vt[i + 1][j + 1];
			}
            
			v[i + 1][jb] = vt[i + 1][jb];
            
			for (int j = jb + 1; j < p_step - 1; j++) {
				v[i][j] = 0.5 * dt * (sigma * sigma * j - (r - q)) * j * v[i + 1][j - 1] + 1 - dt * (sigma * sigma * j * j + r) * v[i + 1][j] + 0.5 * dt * (sigma * sigma * j + (r - q)) * j * v[i + 1][j + 1];
			}

			//coupon
			for (int j = jb; j < p_step - 1; j++) {
				vt[i][j] = vt[i][j] + coupon * coupon_no * 1000 * exp(-r * 3.0 / 365);
				v[i][j] = v[i][j] + coupon * coupon_no * 1000 * exp(-r * 3.0 / 365);
			}
			coupon_no++;
		}

		// review dates
		if (i == review_dates[1] || i == review_dates[2] || i == review_dates[3]) {
			for (int j = 1; j < j0; j++) {
				vt[i][j] = 0.5 * dt * (sigma * sigma * j - (r - q)) * j * vt[i + 1][j - 1] + 1 - dt * (sigma * sigma * j * j + r) * vt[i + 1][j] + 0.5 * dt * (sigma * sigma * j + (r - q)) * j * vt[i + 1][j + 1];
			}

			for (int j = j0; j < p_step - 1; j++) {
				vt[i][j] = 1000 * exp(-r * 3.0 / 365);
			}

			v[i + 1][jb] = vt[i + 1][jb];

			for (int j = jb + 1; j < j0; j++) {
				v[i][j] = 0.5 * dt * (sigma * sigma * j - (r - q)) * j * v[i + 1][j - 1] + 1 - dt * (sigma * sigma * j * j + r) * v[i + 1][j] + 0.5 * dt * (sigma * sigma * j + (r - q)) * j * v[i + 1][j + 1];
			}

			for (int j = j0; j < p_step - 1; j++) {
				v[i][j] = 1000 * exp(-r * 3.0 / 365);
			}

			for (int j = jb; j < p_step - 1; j++) {
				vt[i][j] = vt[i][j] + coupon * coupon_no * 1000 * exp(-r * 3.0 / 365);
				v[i][j] = v[i][j] + coupon * coupon_no * 1000 * exp(-r * 3.0 / 365);
			}
			coupon_no++;
		}

		//regular
		for (int j = 1; j < p_step - 1; j++) {
			vt[i][j] = a = 0.5 * dt * (sigma * sigma * j - (r - q)) * j * vt[i + 1][j - 1] + 1 - dt * (sigma * sigma * j * j + r) * vt[i + 1][j] + 0.5 * dt * (sigma * sigma * j + (r - q)) * j * vt[i + 1][j + 1];
		}

		v[i + 1][jb] = vt[i + 1][jb];

		for (int j = jb + 1; j < p_step - 1; j++) {
			v[i][j] = a = 0.5 * dt * (sigma * sigma * j - (r - q)) * j * v[i + 1][j - 1] + 1 - dt * (sigma * sigma * j * j + r) * v[i + 1][j] + 0.5 * dt * (sigma * sigma * j + (r - q)) * j * v[i + 1][j + 1];
		}
	}
	cout << v[0][j0 + 1] << endl;
}
