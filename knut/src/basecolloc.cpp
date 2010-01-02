// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#include <cmath>
#include "basecolloc.h"
#include "polynomial.h"
#include "system.h"

#ifdef DEBUG
#include <iostream>
#include <fstream>
#endif /*DEBUG*/

static void equidist(Vector& mesh)
{
  const int m = mesh.size();
  for (int i = 0; i < m; i++) mesh(i) = (double)i / ((double)m - 1);
}

static const double legendre_roots[25][25] = 
  {{0},
  {0.5}, // 1
  {0.21132486540518713,0.7886751345948129}, // 2
  {0.1127016653792583,0.5,0.8872983346207417}, // 3
  {0.06943184420297371,0.33000947820757187,0.6699905217924281,0.9305681557970262}, // 4
  {0.04691007703066796,0.2307653449471585,0.5,0.7692346550528415,0.953089922969332}, // 5
  {0.03376524289842403,0.1693953067668677,0.3806904069584015,0.6193095930415985,
    0.8306046932331324,0.9662347571015761}, // 6
  {0.025446043828620757,0.12923440720030277,0.2970774243113014,0.5,
    0.7029225756886986,0.8707655927996972,0.9745539561713793}, // 7
  {0.01985507175123158,0.10166676129318652,0.23723379504183545,0.4082826787521751,
    0.591717321247825,0.7627662049581645,0.898333238706813,0.980144928248768}, // 8
  {0.015919880246187068,0.0819844463366825,0.19331428364970366,0.3378732882980955,
    0.5,0.6621267117019045,0.806685716350295,0.9180155536633178,0.9840801197538123}, // 9
  {0.013046735741415516,0.06746831665550812,0.16029521585048856,0.2833023029353763,
    0.4255628305091844,0.5744371694908156,0.7166976970646237,0.8397047841495116,
    0.9325316833444912,0.9869532642585845}, // 10
  { 0.010885670926971014,  0.056468700115950565,0.13492399721297416,0.24045193539659415,
    0.3652284220238275,0.5,0.6347715779761725,0.7595480646034056,0.8650760027870263,
    0.9435312998840468,0.9891143290730284}, // 11
  { 0.009219682876644486, 0.04794137181476521, 0.11504866290284893, 0.20634102285669104,
    0.3160842505009099,   0.43738329574426554, 0.5626167042557344,  0.6839157494990902,
    0.7936589771433094,   0.8849513370971523,  0.9520586281852412,  0.9907803171233556 }, // 12
  { 0.007908472640693165, 0.041200800388523084,0.09921095463334423, 0.17882533027983172,
    0.2757536244817764,   0.3847708420224326,  0.5,                 0.6152291579775674,
    0.7242463755182234,   0.8211746697201705,  0.900789045366654,   0.9587991996114867,
    0.9920915273592984}, // 13
  {0.006858095651593732,  0.03578255816821063, 0.08639934246513187,0.15635354759415493,
    0.24237568182092356,0.34044381553605507,0.44597252564632817,0.5540274743536718,
    0.6595561844639449,0.757624318179076,0.8436464524058471,0.9136006575348972,
    0.9642174418317947,0.993141904348404},
  {0.0060037409897608085,0.031363303799564035,0.07589670829481437,0.1377911343199068,
    0.21451391369571127,0.30292432646121836,0.39940295300128276,0.5,0.6005970469987173,
    0.6970756735387815,0.7854860863042676,0.862208865680093,0.9241032917051838,
    0.9686366962003857,0.9939962590102489},
  {0.005299532504252191,0.027712488463208285,0.0671843988061569,0.12229779582246603,
    0.19106187779868644,0.2709916111713858,0.35919822461037054,0.4524937450811813,
    0.5475062549188188,0.6408017753896295,0.7290083888286145,0.8089381222013095,
    0.8777022041775406,0.932815601193808,0.9722875115365677,0.9947004674958182},
  {0.004712262342856932,0.024662239115523232,0.059880423136494776,0.10924299805159132,
    0.17116442039165874,0.24365473145676136,0.3243841182730616,0.41075790925207606,
    0.5,0.5892420907479239,0.6756158817269382,0.7563452685432404,0.8288355796083402,
    0.8907570019484374,0.9401195768633666,0.9753377608845084,0.9952877376571484},
  {0.004217415789742551,0.022088025214056728,0.053698766756720695,0.09814752051370701,
    0.15415647846984937,0.2201145844630229,0.29412441926857913,0.37405688715424706,
    0.45761249347913235,0.5423875065208676,0.6259431128457528,0.7058755807314211,
    0.7798854155369719,0.8458435215301504,0.9018524794863332,0.9463012332487712,
    0.9779119747865042,0.995782584210471},
  {0.003796578078055335,0.01989592393262496,0.048422048191554656,0.08864267173118828,
    0.13951691133228333,0.19972734766916367,0.26771462931201756,0.34171795001818517,
    0.4198206771798873,0.5,0.5801793228201126,0.6582820499818149,0.7322853706879799,
    0.8002726523307858,0.8604830886676103,0.9113573282684008,0.951577951807693,
    0.9801040760671133,0.9962034219221172},
  {0.003435700409276654,0.018014036359358776,0.043882785874309604,0.0804415140882958,
    0.12683404677031385,0.18197315963683747,0.24456649902458294,0.3131469556422919,
    0.3861070744291774,0.46173673943325133,0.5382632605667487,0.6138929255708225,
    0.6868530443577096,0.7554335009754155,0.8180268403630893,0.8731659532303624,
    0.9195584859099291,0.9561172141249067,0.9819859636417109,0.996564299593517},
  {0.0031239146916462457,0.01638658071239918,0.039950332934399646,0.07331831770934533,
    0.11578001826292433,0.1664305979012391,0.2241905820563811,0.28782893989628267,
    0.355989341598799,0.42721907291955247,0.5,0.5727809270804476,0.6440106584012006,
    0.7121710601037169,0.7758094179436078,0.8335694020986614,0.884219981737923,
    0.9266816822937112,0.9600496670770324,0.9836134192858432,0.9968760853113308},
  {0.002852707248337527,0.01496975111014015,0.036521613894613225,0.06709371114172646,
    0.10609159700944465,0.15275636840682716,0.20617979824657523,0.2653220810066297,
    0.3290320895539578,0.3960697866558893,0.4651303633401389,0.5348696366598611,
    0.6039302133441106,0.6709679104460418,0.7346779189933691,0.7938202017535103,
    0.8472436315936023,0.893908402989533,0.9329062888555137,0.963478386108585,
    0.9850302489049136,0.9971472927476741},
  {0.002615332496858358,0.013728764391654769,0.033514456591435926,0.06162382086690804,
    0.09755579919043539,0.14066931843746722,0.19019506211810377,0.24524926107696787,
    0.3048494809848553,0.3679321595148276,0.43337158785076696,0.5,0.5666284121492331,
    0.6320678404851723,0.6951505190151435,0.7547507389230321,0.8098049378817609,
    0.8593306815708692,0.9024442008157931,0.9383761791245382,0.9664855433944164,
    0.9862712356102815,0.997384667495246},
  {0.0024063900424206586,0.012635722008147443,0.030862724035082723,0.05679223650356424,
    0.08999900700735042,0.1299379042070915,0.1759531740314857,0.22728926430555285,
    0.2831032461869818,0.3424786601519182,0.40444056626319186,0.4679715535686972,
    0.5320284464313028,0.5955594337368082,0.6575213398480825,0.7168967538130162,
    0.772710735694192,0.824046825968634,0.8700620957880745,0.9100009930382996,
    0.9432077635144789,0.9691372760420776,0.9873642779936926,0.9975936099955564}};

static void poly_gau(Vector& roots)
{
  const int m = roots.size();
  /* construct the matrix */
  P_ERROR_X3(m < static_cast<int>(sizeof(legendre_roots[0])/sizeof(double)), "Unavailable polynomial order ", m, ".");
  for (int i = 0; i < m; ++i)
  {
    roots(i) = legendre_roots[m][i];
  }
}

inline static void repr_mesh(Vector& V)
{
  equidist(V);
}

inline static void col_mesh(Vector& V)
{
  poly_gau(V);
}

static inline void poly_mul(Vector& pp, double bb, double aa)
{
  P_ASSERT_X1(pp(pp.size() - 1) == 0.0, "poly_linmul: truncating the highest order term!");
  for (int i = pp.size() - 1; i > 0; --i)
  {
    pp(i) = aa * pp(i) + bb * pp(i - 1);
  }
  pp(0) = aa * pp(0);
}

static void poly_diff_int_trap(Matrix& sum, const Vector& t)
{
  const int dim = t.size();
  const int DIVS = 8*dim;
  sum.Clear();
  Vector lout(dim);
  Vector dout(dim);
  double d;
  for (int k = 0; k < DIVS; ++k)
  {
    poly_lgr(t, lout, k*1.0/((double)DIVS-1.0));
    poly_dlg(t, dout, k*1.0/((double)DIVS-1.0));
    if ((k==0)||(k==DIVS-1)) d = 0.5/(DIVS-1);
    else d = 1.0/(DIVS-1);
    for (int i = 0; i < dim; i++)
	{
	  for (int j = 0; j < dim; j++)
	  {
	    sum(i,j) += d*dout(i)*lout(j);
	  }
	}
  }
}

static void poly_int_trap(Matrix& sum, const Vector& t)
{
  const int dim = t.size();
  const int DIVS = 8*dim;
  sum.Clear();
  Vector lout(dim);
  double d;
  for (int k = 0; k < DIVS; ++k)
  {
    poly_lgr(t, lout, k*1.0/((double)DIVS-1.0));
    if ((k==0)||(k==DIVS-1)) d = 0.5/(DIVS-1);
    else d = 1.0/(DIVS-1);
    for (int i = 0; i < dim; i++)
	{
	  for (int j = 0; j < dim; j++)
	  {
	    sum(i,j) += d*lout(i)*lout(j);
	  }
	}
  }
}

static void poly_int(Matrix& out, const Vector& t)
{
  int i, j, k;
  Vector poly(2*t.size());
  Vector poly1(t.size());
  Vector poly2(t.size());

  for (i = 0; i < t.size(); i++)
  {
    for (j = 0; j < t.size(); j++)
    {
      poly.Clear();
      poly1.Clear();
      poly2.Clear();
      poly1(0) = 1.0;
      poly2(0) = 1.0;
      //      poly.Print();
      // i,j az out matrix indexe
      //      cout<<"in:poly_mul\n";
      for (k = 0; k < t.size(); k++)
      {
        if (k != i)
          poly_mul(poly1, 1.0, -t(k));
        if (k != j)
          poly_mul(poly2, 1.0, -t(k));
      }
      poly1 /= poly_eval(poly1,t(i));
      poly2 /= poly_eval(poly2,t(j));
      poly_coeff_mul(poly, poly1, poly2);
      // integrate
      for (k = 0; k < poly.size(); k++) poly(k) /= k + 1.0;
      out(i, j) = 0.0;
      // evaluate at x = 0..1
      for (k = 0; k < poly.size(); k++) out(i, j) += poly(k);
    }
  }
//   out.Clear();
//   for (i = 0; i < t.size(); i++) out(i,i) = 1.0/t.size();
}

static void poly_diff_int(Matrix& out, const Vector& t)
{
  int i, j, k;
  Vector poly(2*t.size());
  Vector poly1(t.size());
  Vector poly1_d(t.size());
  Vector poly2(t.size());

  for (i = 0; i < t.size(); i++)
  {
    for (j = 0; j < t.size(); j++)
    {
      poly.Clear();
      poly1.Clear();
      poly1_d.Clear();
      poly2.Clear();
      poly1(0) = 1.0;
      poly2(0) = 1.0;
      //      poly.Print();
      // i,j az out matrix indexe
      //      cout<<"in:poly_mul\n";
      for (k = 0; k < t.size(); k++)
      {
        if (k != i)
          poly_mul(poly1, 1.0, -t(k));
        if (k != j)
          poly_mul(poly2, 1.0, -t(k));
      }
      poly1 /= poly_eval(poly1,t(i));
      poly2 /= poly_eval(poly2,t(j));
      poly_coeff_diff(poly1_d, poly1);
      poly_coeff_mul(poly, poly1_d, poly2);
      // integrate
      for (k = 0; k < poly.size(); k++) poly(k) /= k + 1.0;
      out(i, j) = 0.0;
      // evaluate at x = 0..1
      for (k = 0; k < poly.size(); k++) out(i, j) += poly(k);
    }
  }
}

#define NDIM ndim
#define NPAR npar
#define NINT nint
#define NDEG ndeg

PerSolColloc::PerSolColloc(System& _sys, const int _nint, const int _ndeg) :
    ndim(_sys.ndim()), npar(_sys.npar()),
    nint(_nint), ndeg(_ndeg),
    time(nint*ndeg),
    timeMSH(ndeg*nint + 1),
    metric(ndeg + 1, ndeg + 1),
    metricPhase(ndeg + 1, ndeg + 1),
    mesh(nint + 1), 
    meshINT(ndeg + 1),
    col(ndeg),
    lgr(ndeg+1, ndeg+1),
    out(ndeg + 1)
{
  sys = &_sys;
  for (int i = 0; i < nint + 1; i++) mesh(i) = i * 1.0 / nint;
  repr_mesh(meshINT);
  col_mesh(col);
  
  poly_int_trap(metric, meshINT);
  poly_diff_int_trap(metricPhase, meshINT);
//  meshINT.Print();
//  col.Print();
//  metric.Print();
//  metricPhase.Print();

  // computes the largrange coefficients
  for (int i = 0; i < ndeg+1; i++)
  {
    poly_coeff_lgr(lgr(i), meshINT, i);
  }
}

int PerSolColloc::meshlookup(const Vector& mesh, double t)
{
  // binary search for in which interval is t-tau(k)
  int mid, low = 0, up = mesh.size() - 1;
  while (up - low > 1)
  {
    mid = low + (up - low) / 2;
    if ((mesh(low) <= t) && (mesh(mid) > t)) up = mid;
    else low = mid;
  }
  return low;
}

//----------------------------------------------------------------------------
//
// here are the mesh adaptation routines
//
//----------------------------------------------------------------------------

static void meshConstruct(Vector& newmesh, const Vector& oldmesh, const Vector& eqf)
{
//   for (int i = 1; i < eqf.size()-1; i++) if (isnan(eqf(i))) std::cout<<i<<": nan ";
//   std::cout<<"first "<<eqf(1)<<" end "<<eqf(NINT)<<" ratio "<< eqf(1)/eqf(NINT)<<"\n";
  // now computing the new mesh
  const int nint = oldmesh.size()-1;
  newmesh(0) = 0.0;
  for (int i = 1; i < newmesh.size()-1; i++)
  {
    const double t = eqf(nint)*i/(newmesh.size()-1);
    const int idx = PerSolColloc::meshlookup( eqf, t );
    const double d = (t - eqf(idx))/(eqf(idx+1)-eqf(idx));
//     std::cout<<t<<":"<<d<<":"<<i<<":"<<idx<<":"<<mesh(idx) + d*(mesh(idx+1)-mesh(idx))<<" : "<<mesh(idx+1)-mesh(idx)<<"\n";
    newmesh(i) = oldmesh(idx) + d*(oldmesh(idx+1)-oldmesh(idx));
//     if (eqf(i) < eqf(i-1)) std::cout<<"bad "<<eqf(i-1)<<", "<<eqf(i)<<"\n";
//     if (newmesh(i) < newmesh(i-1)) std::cout<<"very bad "<<newmesh(i-1)<<", "<<newmesh(i)<<"\n";
  }
  newmesh(newmesh.size()-1) = 1.0;
}

static void meshAssess(Vector& eqf, const Vector& mesh, const Vector& profile, const Array1D< Array1D<double> >& lgr)
{
  int ndeg_ = lgr.size() - 1;
  int nint_ = mesh.size() - 1;
  P_ERROR_X1( profile.size() % (ndeg_*nint_ + 1) == 0, "Wrong profile size.");
  int ndim_ = profile.size() / (ndeg_*nint_ + 1);
  
  // compute the coeff of the highest degree term in each interval
  bool small_deri = true;
  const double hmach = 1e-6;
  Matrix hd(nint_+1, ndim_);
  for (int i = 0; i < nint_; i++)
  {
    for (int p = 0; p < ndim_; p++)
    {
      hd(i, p) = 0.0;
      for (int j = 0; j < ndeg_+1; j++)
      {
        hd(i,p) += lgr(j)(ndeg_)*profile( p + ndim_*( j + ndeg_*i ) );
      }
      // adjust by the mesh interval
      hd(i,p) /= pow(mesh(i+1)-mesh(i),ndeg_);
      if (fabs(hd(i,p)) > hmach) small_deri = false;
    }
  }
//   if (small_deri) std::cout<<"small derivatives\n";
  // takes care of periodicity
  // this has to be changed when other boundary condition is used
  for (int p = 0; p < ndim_; p++)
  {
    hd(nint_,p) = hd(0,p);
  }
  // computes the (m+1)-th derivative.
  // The mesh modulo need not to be changed, when not periodic BC is used.
  for (int i = 0; i < nint_; i++)
  {
    double dtav;
    if ( i+2 < nint_ ) dtav = 0.5*(mesh(i+2)-mesh(i));
    else dtav = 0.5*(1.0+mesh((i+2)-nint_)-mesh(i));
    if( dtav < 0.0 ) std::cout<<"dtav<0\n";
    for (int p = 0; p < ndim_; p++)
    {
      hd(i,p) = (hd(i+1,p) - hd(i,p))/dtav;
    }
  }
  // takes care of periodicity
  // this has to be changed when other boundary condition is used
  for (int p = 0; p < ndim_; p++)
  {
    hd(nint_,p) = hd(0,p);
  }
  // eqf contains the integral which has to be equidistributed
  P_ERROR_X1( eqf.size() == nint_+1, "EQF has wrong size.");
  // when the derivatives are too small;
  eqf(0) = 0.0;
  // computing eqf
  const double pwr=1.0/(ndeg_+1.0);
  for (int j=0; j < nint_; ++j)
  {
    double EP=0;
    if (j == 0)
    {
      for (int i = 0; i < ndim_; ++i)
      {
        EP+=pow(fabs(hd(nint_,i)),pwr);
      }
    } else
    {
      for (int i = 0; i < ndim_; ++i)
      {
        EP+=pow(fabs(hd(j-1,i)),pwr);
      }
    }
    double E=0;
    for (int i = 0; i < ndim_; ++i)
    {
      E+=pow(fabs(hd(j,i)),pwr);
    }
    eqf(j+1)=eqf(j)+0.5*(mesh(j+1)-mesh(j))*(E+EP);
    if (E+EP==0) small_deri = true;
  }
  if (small_deri) for (int i = 0; i < nint_+1; ++i) eqf(i) = i;
}

void PerSolColloc::meshAdapt_internal( Vector& newmesh, const Vector& profile )
{
  Vector eqf(NINT+1);
  meshAssess(eqf, mesh, profile, lgr);
  meshConstruct(newmesh, mesh, eqf);
}

static void profileConvert(Vector& newprofile, const Vector& newmesh, const Vector& profile, const Vector& mesh,
                           const Array1D< Array1D<double> >& old_lgr, const int ndim)
{
  const int old_nint = mesh.size()-1;
  const int old_ndeg = (profile.size()/ndim - 1)/old_nint;
  const int new_nint = newmesh.size()-1;
  const int new_ndeg = (newprofile.size()/ndim - 1)/new_nint;
  // creating the new profile with the same old meshINT
  for (int i = 0; i < new_nint; i++)
  {
    for (int j = 0; j < new_ndeg; j++)
    {
      const double t = newmesh(i) + j*(newmesh(i+1)-newmesh(i))/new_ndeg;
      int idx = PerSolColloc::meshlookup( mesh, t );
      const double d = (t - mesh(idx))/(mesh(idx+1)-mesh(idx));
      for (int p = 0; p < ndim; p++)
        newprofile(p+ndim*(j+i*new_ndeg)) = 0.0;
      for ( int k = 0; k < old_ndeg+1; k++)
      {
        double c = poly_eval(old_lgr(k), d);
        for (int p = 0; p < ndim; p++)
        {
          newprofile(p+ndim*(j+i*new_ndeg)) += c*profile(p+ndim*(k+idx*old_ndeg));
        }
      }
    }
  }
  for (int p = 0; p < ndim; p++)
    newprofile(p+ndim*(new_ndeg*(newmesh.size()-1))) = profile(p+ndim*(old_ndeg*(mesh.size()-1)));
}

void PerSolColloc::meshAdapt(Vector& newprofile, const Vector& profile, Vector& newtangent, const Vector& tangent)
{
  // saving the solution into solNu
#ifdef DEBUG
  Vector profile_tmp(profile);
#endif

  Vector newmesh(mesh);
  meshAdapt_internal(newmesh, profile);
  profileConvert(newprofile, newmesh, profile, mesh, lgr, NDIM);
  profileConvert(newtangent, newmesh, tangent, mesh, lgr, NDIM);
#ifdef DEBUG
  // printing the adapted profile
  std::ofstream file1("prof");
  file1 << std::scientific;
  file1.precision(12);
  std::ofstream file2("newprof");
  file2 << std::scientific;
  file2.precision(12);

  std::ofstream file3("gradient");
  file3 << std::scientific;
  file3.precision(12);
  for (int i=0; i<NINT; ++i)
  {
    for (int j=0; j<NDEG+1; ++j)
    {
      const double t1 = mesh(i) + j*(mesh(i+1)-mesh(i))/NDEG;
      const double t2 = newmesh(i) + j*(newmesh(i+1)-newmesh(i))/NDEG;
      file1<<t1<<"\t";
      file2<<t2<<"\t";
      file3<<t1<<"\t"<<t2<<"\t";
      for (int p=0; p<NDIM; ++p)
      {
	file1<<profile_tmp(p+NDIM*(j+NDEG*i))<<"\t";
	file2<<profile(p+NDIM*(j+NDEG*i))<<"\t";
      }
      file1<<"\n";
      file2<<"\n";
      file3<<"\n";
    }
  }
#endif //DEBUG
  mesh = newmesh;
}

//----------------------------------------------------------------------------
//
// here are the integration routines
//
//----------------------------------------------------------------------------

void PerSolColloc::getMetric(Matrix& mt, const Vector& t)
{
  poly_int(mt, t);
}

void PerSolColloc::getDiffMetric(Matrix& mt, const Vector& t)
{
  poly_int(mt, t);
}

void PerSolColloc::star(Vector& out, const Vector& in, const Matrix& mt, const Vector& msh, int dim)
{
  const int t_deg = mt.Col() - 1;
  const int t_int = (out.size() / dim - 1) / t_deg;
  P_ERROR(in.size() == out.size());
  P_ERROR(out.size() == dim*(t_deg*t_int + 1));
  P_ERROR(msh.size() == t_int + 1);
  out.Clear();
  for (int i = 0; i < t_int; ++i)
  {
    const double dx = msh(i + 1) - msh(i);
    // ez itt a matrixszorzas
    for (int k = 0; k < t_deg + 1; ++k)
    {
      for (int l = 0; l < t_deg + 1; ++l)
      {
        for (int j = 0; j < dim; ++j)
        {
          out(j + dim*(k + i*t_deg)) += dx * mt(k, l) * in(j + dim * (l + i * t_deg));
        }
      }
    }
  }
#ifdef MADD // whether we need to add the headpoint ?
  for (int j = 0; j < dim; j++)
  {
    out(t_int*t_deg*dim + j) += in(t_int * t_deg * dim + j);
  }
#endif
}

double PerSolColloc::integrate(const Vector& v1, const Vector& v2, const Matrix& mt, const Vector& msh, int dim)
{
  double res = 0.0;
  const int t_deg = mt.Col() - 1;
  const int t_int = (v1.size() / dim - 1) / t_deg;
  P_ERROR(v1.size() == v2.size());
  P_ERROR(v1.size() == dim*(t_deg*t_int + 1));
  P_ERROR(msh.size() == t_int + 1);
  for (int i = 0; i < t_int; ++i)
  {
    const double dx = msh(i + 1) - msh(i);
    // ez itt a matrixszorzas
    for (int k = 0; k < t_deg; ++k)
    {
      for (int l = 0; l < t_deg; ++l)
      {
        for (int j = 0; j < dim; ++j)
        {
          res += v1(j + dim * (k + i * t_deg)) * dx * mt(k, l) * v2(j + dim * (l + i * t_deg));
        }
      }
    }
  }
#ifdef MADD // whether we need to add the headpoint ?
  for (int j = 0; j < dim; j++)
  {
    res += v1(t_int * t_deg * dim + j) * v2(t_int * t_deg * dim + j);
  }
#endif
  return res;
}

void PerSolColloc::Star(Vector& V1, const Vector& V2)
{
  V1.Clear();
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          V1(j + NDIM*(k + i*NDEG)) += dx * metric(k, l) * V2(j + NDIM * (l + i * NDEG));
        }
      }
    }
  }

#ifdef MADD // whether we need to add the headpoint ?

  // Now we add (M udot)^* M\phi + <udot, \phi > ...
  for (int j = 0; j < NDIM; j++)
  {
    V1(NINT*NDEG*NDIM + j) += V2(NINT * NDEG * NDIM + j);
  }

#endif
}

double PerSolColloc::Integrate(const Vector& V1, const Vector& V2)
{
  double res = 0.0, head = 0.0;
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          res += dx * V1(j + NDIM * (k + i * NDEG)) * metric(k, l) * V2(j + NDIM * (l + i * NDEG));
        }
      }
    }
  }

#ifdef MADD
  for (int j = 0; j < NDIM; j++)
  {
    head += V1(NINT * NDEG * NDIM + j) * V2(NINT * NDEG * NDIM + j);
  }
#endif

  return res + head;
}

double PerSolColloc::IntegrateCont(const Vector& V1, const Vector& V2, const Vector& V3)
{
  double res = 0.0, head = 0.0;
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          res += dx * V1(j + NDIM * (k + i * NDEG)) * metric(k, l) * (V2(j + NDIM * (l + i * NDEG)) - V3(j + NDIM * (l + i * NDEG)));
        }
      }
    }
  }
#ifdef MADD
  for (int j = 0; j < NDIM; j++)
  {
    head += V1(NINT * NDEG * NDIM + j) * (V2(NINT * NDEG * NDIM + j) - V3(NINT * NDEG * NDIM + j));
  }
#endif
  return res + head;
}

void PerSolColloc::PhaseStar(Vector& V1, const Vector& V2)
{
  V1.Clear();
  for (int i = 0; i < NINT; i++)
  {
    for (int j = 0; j < NDIM; j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          V1(j + NDIM*(k + i*NDEG)) += metricPhase(k, l) * V2(j + NDIM * (l + i * NDEG));
        }
      }
    }
  }
}

void PerSolColloc::PhaseRotStar(Vector& V1, const Vector& V2, const Array1D<int>& Re, const Array1D<int>& Im)
{
  V1.Clear();
  for (int i = 0; i < NINT; i++)
  {
    const double dx = mesh(i + 1) - mesh(i);
    for (int j = 0; j < Re.size(); j++)
    {
      // ez itt a matrixszorzas
      for (int k = 0; k < NDEG + 1; k++)
      {
        for (int l = 0; l < NDEG + 1; l++)
        {
          V1(Re(j) + NDIM*(k + i*NDEG)) -= dx * metric(k, l) * V2(Im(j) + NDIM * (l + i * NDEG));
          V1(Im(j) + NDIM*(k + i*NDEG)) += dx * metric(k, l) * V2(Re(j) + NDIM * (l + i * NDEG));
        }
      }
    }
  }
}

void PerSolColloc::pdMeshConvert(Vector& newprofile, Vector& newtangent, const Vector& oldprofile, const Vector& oldtangent)
{
  Vector tmp_mesh(2*NINT+1);
  Vector tmp_profile(NDIM*(2*NINT*NDEG+1));
  Vector tmp_tangent(NDIM*(2*NINT*NDEG+1));
  for (int i = 0; i < NINT; ++i)
  {
    tmp_mesh(i) = 0.0 + 0.5*mesh(i);
    tmp_mesh(NINT+i) = 0.5 + 0.5*mesh(i);
    for (int j = 0; j < NDEG; ++j)
    {
      for (int p = 0; p < NDIM; ++p)
      {
        tmp_profile(p+NDIM*(j+i*NDEG))        = oldprofile(p+NDIM*(j+i*NDEG));
        tmp_profile(p+NDIM*(j+(i+NINT)*NDEG)) = oldprofile(p+NDIM*(j+i*NDEG));
        tmp_tangent(p+NDIM*(j+i*NDEG))        = oldtangent(p+NDIM*(j+i*NDEG));
        tmp_tangent(p+NDIM*(j+(i+NINT)*NDEG)) = -oldtangent(p+NDIM*(j+i*NDEG));
      }
    }
  }
  tmp_mesh(2*NINT) = 1.0;
  for (int p = 0; p < NDIM; ++p)
  {
    tmp_profile(p+NDIM*(2*NINT*NDEG)) = tmp_profile(p);
    tmp_tangent(p+NDIM*(2*NINT*NDEG)) = tmp_tangent(p);
  }
  // constructing the new mesh
  Vector eqf(tmp_mesh.size());
  for (int i = 0; i < eqf.size(); ++i) eqf(i) = i;
  meshConstruct(mesh, tmp_mesh, eqf);
  profileConvert(newprofile, mesh, tmp_profile, tmp_mesh, lgr, NDIM);
  profileConvert(newtangent, mesh, tmp_tangent, tmp_mesh, lgr, NDIM);
}

void PerSolColloc::Import(Vector& newprofile, const Vector& oldprofile, const Vector& oldmesh, int old_ndeg, bool adapt)
{
  Vector old_meshINT(old_ndeg+1);
  repr_mesh(old_meshINT);

  Array1D< Array1D<double> > old_lgr(old_ndeg+1, old_ndeg+1);
  for (int i = 0; i < old_ndeg+1; i++)
  {
    poly_coeff_lgr(old_lgr(i), old_meshINT, i);
  }

  Vector eqf(oldmesh.size());
  if (adapt) meshAssess(eqf, oldmesh, oldprofile, old_lgr);
  else for (int i = 0; i < oldmesh.size(); ++i) eqf(i) = i;
  meshConstruct(mesh, oldmesh, eqf);
  profileConvert(newprofile, mesh, oldprofile, oldmesh, old_lgr, NDIM);
}

// it exports for CollocTR and PointTR, so no last value is necessary
void PerSolColloc::Export(Vector& outs, const Vector& mshint, const Vector& mshdeg, const Vector& in)
{
  int nint_ = mshint.size() - 1;
  int ndeg_ = mshdeg.size() - 1;
  Vector in_mesh(NDEG + 1);
  Vector in_lgr(NDEG + 1);

  for (int i = 0; i < NDEG + 1; i++) in_mesh(i) = i * 1.0 / NDEG;

  for (int i = 0; i < nint_; i++)
  {
    for (int j = 0; j < ndeg_; j++)
    {
      double t = mshint(i) + mshdeg(j) / nint_;
      int k = meshlookup(mesh, t);
      // std::cout<<"int "<<i<<" "<<k<<"\n";
      double c = (t - mesh(k)) / (mesh(k + 1) - mesh(k));  // mesh is the interval mesh in the class

      poly_lgr(in_mesh, in_lgr, c);
      // in_lgr.Print();
      for (int p = 0; p < NDIM; p++)
      {
        outs(p + NDIM*(j + i*ndeg_)) = 0.0;
        for (int r = 0; r < NDEG + 1; r++)
        {
          outs(p + NDIM*(j + i*ndeg_)) += in(p + NDIM * (r + k * NDEG)) * in_lgr(r);
        }
      }
    }
  }
}
