
class contParams
{
  public:
    contParams()
        : sysname(),
        label(0),
        pttype(SolTF),
        cp(VarNone),
        parx(0),
        eqns(0),
        vars(0),
        npar(0),
        nint(0),
        ndeg(0),
        nmul(0),
        stab(false),
        nmat(0),
        nint1(0),
        nint2(0),
        ndeg1(0),
        ndeg2(0),
        steps(0),
        cpMin(0),
        cpMax(0),
        ds(0.0),
        dsMin(0.0),
        dsMax(0.0),
        dsStart(0.0),
        epsC(0.0),
        epsR(0.0),
        epsS(0.0),
        nitC(0.0),
        nitR(0.0),
        nitS(0.0),
        symre(0),
        symim(0),
        parxType(0, false),
        eqnsType(0),
        varsType(0, true)
    { }
    void setSysName();
  private:
    std::string sysname;
    int    label;
    PtType pttype;
    Var    cp;
    std::vector<Var> parx;
    std::vector<Eqn> eqns;
    std::vector<Var> vars;
    int    npar;
    int    nint;
    int    ndeg;
    int    nmul;
    bool   stab;
    int    nmat;
    int    nint1;
    int    nint2;
    int    ndeg1;
    int    ndeg2;
    int    steps;
    double cpMin;
    double cpMax;
    double ds;
    double dsMin;
    double dsMax;
    double dsStart;
    double epsC;
    double epsR;
    double epsS;
    int    nitC;
    int    nitR;
    int    nitS;
    std::vector<int> symre;
    std::vector<int> symim;
    VarType      parxType;
    EqnType      eqnsType;
    VarType      varsType;
};
