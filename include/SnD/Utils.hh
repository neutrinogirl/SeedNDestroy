//
// Created by Stephane Zsoldos on 6/27/22.
//

#ifndef WRATTER_INCLUDE_WRATTER_UTILS_HH_
#define WRATTER_INCLUDE_WRATTER_UTILS_HH_

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
 private:
  double RINDEX() const {return mRIndex;}
  double mRIndex = 1.33;
  double SOL_VACUUM() const {return mSoL_vacuum;}
  double mSoL_vacuum = 299.792; // mm/ns;
  double SOL() const {return mSoL;}
  double mSoL = mSoL_vacuum/mRIndex; // mm/ns;
  Csts() = default;
};

#endif //WRATTER_INCLUDE_WRATTER_UTILS_HH_
