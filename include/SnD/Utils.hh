//
// Created by Stephane Zsoldos on 6/27/22.
//

#ifndef SND_INCLUDE_SND_UTILS_HH_
#define SND_INCLUDE_SND_UTILS_HH_

class Csts {
 public:
  Csts(const Csts&) = delete;
  static Csts& Get(){
	static Csts instance;
	return instance;
  }
  static double GetRIndex()     {return Get().RINDEX();};
  static double GetSoL_vacuum() {return Get().SOL_VACUUM();};
  static double GetSoL()        {return Get().SOL();};
  static double GetSqrt2()      {return Get().SQRT2();};
 private:
  [[nodiscard]] double RINDEX() const {return mRIndex;}
  double mRIndex = 1.33;
  [[nodiscard]] double SOL_VACUUM() const {return mSoL_vacuum;}
  double mSoL_vacuum = 299.792; // mm/ns;
  [[nodiscard]] double SOL() const {return mSoL;}
  double mSoL = mSoL_vacuum/mRIndex; // mm/ns;
  [[nodiscard]] double SQRT2() const {return mSqrt2;}
  double mSqrt2 = 1.41421356237; // sqrt(2);
  Csts() = default;
};

class Gauss {
 public:
  Gauss(const Gauss&) = delete;
  static Gauss& Get(){
	static Gauss instance;
	return instance;
  }
  static double Eval(double x, double mu, double sigma) {
	return Get().Eval_(x, mu, sigma);
  }
 private:
  static double Eval_(double x, double mu, double sigma) {
	return 1 / (sigma*Csts::GetSqrt2()*sqrt(M_PI)) * exp(-0.5*pow((x-mu)/sigma, 2));
  }
  Gauss() = default;
};

class IQR {
 public:
  IQR(const IQR&) = delete;
  static IQR& Get(){
	static IQR instance;
	return instance;
  }
  static double GetIQR(const std::vector<double>& data) {
	return Get().GetIQR_(data);
  }
 private:
  static double GetIQR_(const std::vector<double>& data) {
	std::vector<double> sorted_data = data;
	std::sort(sorted_data.begin(), sorted_data.end());

	int n = static_cast<int>(sorted_data.size());
	int q1_index = n / 4;
	int q3_index = (3 * n) / 4;

	double q1 = sorted_data[q1_index];
	double q3 = sorted_data[q3_index];

	return q3 - q1;
  }
  IQR() = default;
};

#endif //SND_INCLUDE_SND_UTILS_HH_
