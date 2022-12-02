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
  double RINDEX() const {return mRIndex;}
  double mRIndex = 1.33;
  double SOL_VACUUM() const {return mSoL_vacuum;}
  double mSoL_vacuum = 299.792; // mm/ns;
  double SOL() const {return mSoL;}
  double mSoL = mSoL_vacuum/mRIndex; // mm/ns;
  double SQRT2() const {return mSqrt2;}
  double mSqrt2 = 1.41421356237; // sqrt(2);
  Csts() = default;
};

#endif //SND_INCLUDE_SND_UTILS_HH_
