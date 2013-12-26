// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 - 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef EXPR_H
#define EXPR_H

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <functional>
#include <ostream>
#include <iostream>

#define KNUT_STACK_INCREMENT 5

namespace ExpTree
{

enum class TokenType : int
{
  Invalid=0,
  Function,
  UnaryPlus,
  UnaryMinus,
  Plus,
  Minus,
  Times,
  Divide,
  Power,
  LeftParen,
  RightParen,
  Comma,
  Equals,
  Define,
  Semicolon,
  Symbol,
  Number,
  Var,
  Par,
  Expr
};

// the pointer is owned somewhere else
class Value
{
public:
  Value() : data(nullptr), skip(1) { }
  ~Value() { }
  Value (const Value& value) { data = value.data; skip = value.skip; }
  double* data;
  size_t skip;
};

class ConstValue
{
public:
  ConstValue() : data(nullptr), skip(1) { }
  ~ConstValue() { }
  ConstValue (const ConstValue& value) { data = value.data; skip = value.skip; }
  ConstValue (const Value& value) { data = value.data; skip = value.skip; }
  const double* data;
  size_t skip;
};

class ValueStack
{
public:
  ValueStack ();
  ValueStack (double *res, size_t w) : width(w), stack(1)
  {
    stack[0].data = res;
    stack[0].skip = 1;
  }
  ~ValueStack ()
  {
    for (size_t k = 1; k < stack.size(); k++) delete[] stack [k].data;
  }
  void setResult (double *res, size_t skip = 1)
  {
    stack[0].data = res;
    stack[0].skip = skip;
  }
  Value& operator[](size_t p) { return stack[p]; }
  void resizeDepth (size_t depth)
  {
    const size_t from = stack.size();
    if (depth > from) stack.resize(depth);
    for (size_t k = from; k < depth; k++)
    {
      stack [k].data = new double[width];
      stack [k].skip = 1;
    }
  }
  void resizeWidth (size_t w)
  {
    if (w > width)
    {
      for (size_t k = 1; k < stack.size(); k++)
      {
        delete[] stack [k].data;
        stack [k].data = new double[w];
        stack [k].skip = 1;
      }
      width = w;
    }
  }
  void resize (size_t depth, size_t w)
  {
    const size_t from = stack.size();
    const size_t nw = (w > width) ? w : width;
    if (depth > from) stack.resize(depth);
    for (size_t k = from; k < depth; k++)
    {
      stack [k].data = new double[nw];
      stack [k].skip = 1;
    }
    if (w > width)
    {
      for (size_t k = 1; k < from; k++)
      {
        delete[] stack [k].data;
        stack [k].data = new double[w];
        stack [k].skip = 1;
      }
      width = w;
    }
  }
  size_t size() const { return stack.size(); }
private:
  size_t width;
  std::vector<Value> stack;
};

// forward declaration
class NodeSymbol;
class NodeFunction;

class Node
{
public:
  Node (TokenType tp) : type(tp) { }
  virtual ~Node() { }
  virtual void deleteTree () { };
  virtual Node* copy () const = 0;
  virtual size_t children () const = 0;
  virtual const Node* getChild (size_t n) const = 0;
  virtual Node* getChild (size_t n) = 0;
  virtual size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const = 0;
//  virtual void addArgument (size_t pos, Node* arg) = 0;
  virtual void print (std::ostream& out) const = 0;
  virtual bool operator != (const Node& node) const = 0;
  bool operator == (const Node& node) const { return !((*this) != node); }
  virtual bool operator < (const Node& node) const = 0;
  virtual void replaceSymbol(const Node& sym, const Node& node, Node** parent) = 0;
  virtual void replaceSymbol(const Node& node, Node** parent) = 0;
  virtual void findFunction(std::vector<const NodeFunction*>& lst, const std::string& name) const = 0;
  virtual Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const = 0;
  virtual void optimize (Node** parent, const std::function<bool(const Node*)>& zero) = 0;
  virtual void find (TokenType tp, const std::function<void(const Node*)>& found) const;
  void evaluateWithPar (double *res, const std::vector<double>& parInit);
  const TokenType type;
  static const std::function<bool(const Node*)> alwaysFalse;
};

class NodeNumber : public Node
{
public:
  NodeNumber (double val) : Node(TokenType::Number), value(val) {}
  Node* copy () const { return new NodeNumber(value); }
  size_t children () const { return 0; }
  const Node* getChild (size_t n) const { return 0; }
  Node* getChild (size_t n) { return 0; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const { out << value; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) {}
  void replaceSymbol(const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero) {}
  double getValue() const { return value; }
private:
  double value;
};

class NodeSymbol : public Node
{
public:
  NodeSymbol (const std::string& nm) : Node(TokenType::Symbol), name(nm) {}
  Node* copy () const { return new NodeSymbol(name); }
  size_t children () const { return 0; }
  const Node* getChild (size_t n) const { return nullptr; }
  Node* getChild (size_t n) { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const { out << name; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  // replaces the symbol with the node
  // parent points to the pointer of the current class
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol (const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero) {}
  const std::string& getName () const { return name; }
private:
  std::string name;
};

class NodeVar : public Node
{
public:
  NodeVar (size_t id) : Node(TokenType::Var), idx(id) { }
  ~NodeVar() { }
  void deleteTree () {}
  Node* copy () const { return new NodeVar(idx); }
  size_t children () const {return 0; }
  const Node* getChild (size_t n) const { return nullptr; }
  Node* getChild (size_t n) { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const { out << "V["<<idx<<"]"; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol (const Node& sym, const Node& node, Node** parent);
  void replaceSymbol (const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero) {}
  size_t getIdx () const { return idx; }
private:
  size_t idx;
};

class NodePar : public Node
{
public:
  NodePar (size_t id) : Node(TokenType::Par), idx(id) { }
  ~NodePar() { }
  void deleteTree () {}
  Node* copy () const { return new NodePar(idx); }
  size_t children () const { return 0; }
  const Node* getChild (size_t n) const { return nullptr; }
  Node* getChild (size_t n) { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const { out << "P["<<idx<<"]"; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero) {}
  size_t getIdx () const { return idx; }
private:
  size_t idx;
};

class NodeExpr : public Node
{
public:
  NodeExpr (size_t id) : Node(TokenType::Expr), idx(id), np(0), nv(0) { }
  NodeExpr (size_t idx_, size_t np_, size_t nv_, size_t dp_, size_t dv0_, size_t dv1_)
    : Node(TokenType::Expr), idx(idx_), np(np_), nv(nv_), dp(dp_) { dv[0] = dv0_; dv[1] = dv1_; }

  ~NodeExpr() { }
  void deleteTree () {}
  NodeExpr* copy () const
  {
    NodeExpr *ne = new NodeExpr(idx);
    ne->np = np;
    ne->nv = nv;
    ne->dp = dp;
    ne->dv[0] = dv[0];
    ne->dv[1] = dv[1];
    return ne;
  }
  size_t children () const { return 0; }
  const Node* getChild (size_t n) const { return nullptr; }
  Node* getChild (size_t n) { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const
  {
    out << "ex" << idx;
    if (np > 0) out << "_p" << dp;
    if (nv == 1) out << "_v" << dv[0];
    if (nv == 2)
    {
      if (dv[0] < dv[1]) out << "_v" << dv[0] << "_v" << dv[1];
      else out << "_v" << dv[1] << "_v" << dv[0];
    }
  }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero)
  {
    if (zero (this))
    {
      Node* tmp = new NodeNumber (0.0);
      delete *parent;
      *parent = tmp;
      return;
    }
  }
  size_t getIdx () const { return idx; }
  size_t getNp () const { return np; }
  size_t getNv () const { return nv; }
  size_t getDp () const { return dp; }
  size_t getDv0 () const { return dv[0]; }
  size_t getDv1 () const { return dv[1]; }
private:
  size_t idx;
  size_t np;
  size_t nv;
  size_t dp;
  size_t dv[2];
};

class NodeFunction : public Node
{
public:
  NodeFunction (const std::string& nm, size_t nargs) : Node(TokenType::Function), name(nm), args(nargs)
  {
    for (size_t k=1; k<nargs; k++) args[k] = 0;
  }
  void deleteTree() { for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree(); delete args[k]; } }
  Node* copy () const;
  size_t children () const { return args.size(); }
  const Node* getChild (size_t n) const { return args[n]; }
  Node* getChild (size_t n) { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  // replaces the symbol with the node
  // parent points to the pointer of the current class
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero);
  // specific
  void specialize (NodeFunction** parent);
  void addArgument (size_t pos, Node* arg) { args[pos] = arg; }
  const std::string& getName() const { return name; }
protected:
  std::string name;
  std::vector<Node*> args;
};

class NodeAdd : public Node
{
public:
  NodeAdd (size_t nargs) : Node(TokenType::Plus), args(nargs), mul(nargs)
  {
    for (size_t k=1; k<nargs; k++) { args[k] = 0; mul[k] = 1.0; }
  }
  void deleteTree () { for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree (); delete args[k]; } }
  Node* copy () const;
  size_t children () const { return args.size(); }
  const Node* getChild (size_t n) const { return args[n]; }
  Node* getChild (size_t n) { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  // specific
  void addArgument (size_t pos, Node* arg, double m) { args[pos] = arg; mul[pos] = m; }
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero);
  void sort ();
private:
  std::vector<Node*> args;
  std::vector<double> mul;
};

class NodeTimes : public Node
{
public:
  NodeTimes (size_t nargs) : Node(TokenType::Times), args(nargs), smul(1.0), divide(nargs)
  {
    for (size_t k=1; k<nargs; k++) { args[k] = 0; divide[k] = false; }
  }
  void deleteTree () { for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree (); delete args[k]; } }
  Node* copy () const;
  size_t children () const { return args.size(); }
  const Node* getChild (size_t n) const { return args[n]; }
  Node* getChild (size_t n) { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  // needs to handle the divisions
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  // specific
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero);
  void sort ();
  void setMul (double mul) { smul = mul; }
  double getMul () { return smul; }
  void addArgument (size_t pos, Node* arg, double m, bool div);
private:
  std::vector<Node*> args;
  double smul;
  std::vector<bool> divide;
};

class NodePower : public Node
{
public:
  NodePower () : Node(TokenType::Power), args(2) {}
  void deleteTree () { for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree (); delete args[k]; } }
  Node* copy () const;
  size_t children () const { return args.size(); }
  const Node* getChild (size_t n) const { return args[n]; }
  Node* getChild (size_t n) { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const;
  void optimize (Node** parent, const std::function<bool(const Node*)>& zero);
  // specific
  void addArgument (size_t pos, Node* arg) { args[pos] = arg; }
private:
  std::vector<Node*> args;
};

class Expression
{
public:
  Expression () : root (nullptr) { }
  Expression (const Expression& ex) // copy constructor necessary for containers
  {
    if (root != nullptr)
    {
      root -> deleteTree ();
      delete root;
      root = nullptr;
    }
    if (ex.root != nullptr) root = ex.root -> copy ();
  }
  Expression (const std::string& str);
  ~Expression ()
  {
    if (root != nullptr)
    {
      root -> deleteTree ();
      delete root;
    }
  }
  // copy operator
  const Expression& operator = (const Expression& exp); // disable copy operator
  // takes over the management of the Node rt
  void fromString (const std::string& str);
  void fromNode (Node* rt)
  {
    if (root != nullptr)
    {
      root -> deleteTree ();
      delete root;
    }
    root = rt;
  }
  const Node* node () const { return root; }
  void derivative (Expression& ex, const Node* nd,
                   const std::function<bool(const Node*)>& zero = Node::alwaysFalse) const
  {
    if(root != nullptr)
    {
      Node* deri = root->derivative (nd, zero);
      deri -> optimize (&deri, zero);
      ex.fromNode (deri);
    }
  }
  void optimize (const std::function<bool(const Node*)>& zero = Node::alwaysFalse)
  {
    if (root != nullptr) root -> optimize (&root, zero);
  }
  void print (std::ostream& out) const { if (root != nullptr) root -> print (out); }
  Node* copy () const { if(root != nullptr) return root -> copy (); else return nullptr; }
  void evaluate (ValueStack& stack, std::function<void(size_t, const Node*)> fun, size_t sp, size_t len) const
  {
    if (root != nullptr) root -> evaluate (stack, sp, fun, len);
  }
  bool isZero() const
  {
    if (root != nullptr)
    {
      const NodeNumber* num = dynamic_cast<NodeNumber*>(root);
      if (num != nullptr)
      {
        if (std::fabs(num -> getValue()) < 2*std::numeric_limits<double>::epsilon())
        {
          return true;
        } else return false;
      } else return false;
    } else return true;
  }
  size_t test ()
  {
    if (root != nullptr)
    {
      double zero = 0.0;
      double res;
      ValueStack stk (&res, 1);
      auto fun = [&zero, &stk] (size_t sp, const Node*) { stk[sp].data[0] = 0.0; stk[sp].skip = 1; };
      root -> evaluate (stk, 0, fun, 1);
      return stk.size();
    }
    return 0;
  }
  void knutSplit (
    std::string& sysName,
    std::vector<std::string>& varName,
    std::vector<Expression>& varDotExpr,
    std::vector<std::string>& exprName,
    std::vector<Expression>& exprFormula,
    std::vector<Expression>& varInit,
    std::vector<double>& mass,
    std::vector<Expression>& delayExpr,
    std::vector<std::string>& parName,
    std::vector<double>& parInit);
  static bool fromXML (std::string& oexpr, const std::string& xmlstring);
private:
  // takes over ownership
  Expression (Node* rt) : root(rt) { }
  Node* root;
};

} // namespace

#endif
