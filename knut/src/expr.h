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

#include <string>
#include <vector>
#include <functional>
#include <ostream>
#include <iostream>

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
  Semicolon, 
  Symbol, 
  Number, 
  Var,
  Par 
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
  virtual size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const = 0;
  virtual size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const = 0;
//  virtual void addArgument (size_t pos, Node* arg) = 0;
  virtual void print (std::ostream& out) const = 0;
  virtual bool operator != (const Node& node) const = 0;
  bool operator == (const Node& node) const { return !((*this) != node); }
  virtual bool operator < (const Node& node) const = 0;
  virtual void replaceSymbol(const Node& sym, const Node& node, Node** parent) = 0;
  virtual void replaceSymbol(const Node& node, Node** parent) = 0;
  virtual void findFunction(std::vector<const NodeFunction*>& lst, const std::string& name) const = 0;
  virtual Node* derivative (const Node* var) const = 0;
  virtual void optimize (Node** parent) = 0;
  // takes arguments from the stack (as many as necessary
  // after pushing the resulting value on top of the stack
  // returns the resolt on the top of the stack 
  // the return value points to the top of the stack
  // this is for preallocating the stack
  virtual size_t stackCount(size_t sp, size_t& max_stack) = 0;
  const TokenType type;
};

class NodeNumber : public Node
{
public:
  NodeNumber (double val) : Node(TokenType::Number), value(val) {}
  Node* copy () const { return new NodeNumber(value); }
  size_t children () const { return 0; }
  const Node* getChild (size_t n) const { return 0; }
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const { out << value; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) {}
  void replaceSymbol(const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var) const;
  void optimize (Node** parent) {}
  size_t stackCount(size_t sp, size_t& max_stack) { if (max_stack < sp+1) max_stack = sp+1; return sp; }
  double getValue() { return value; }
private:
  double value;
};

class NodeSymbol : public Node
{
public:
  NodeSymbol (const std::string& nm) : Node(TokenType::Symbol), name(nm) {}
  Node* copy () const { return new NodeSymbol(name); }
  size_t children () const { return 0; }
  const Node* getChild (size_t n) const { return 0; }
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const { out << name; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  // replaces the symbol with the node
  // parent points to the pointer of the current class
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol (const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var) const;
  void optimize (Node** parent) {}
  size_t stackCount(size_t sp, size_t& max_stack) { if (max_stack < sp+1) max_stack = sp+1; return sp; }
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
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const { out << "V["<<idx<<"]"; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol (const Node& sym, const Node& node, Node** parent);
  void replaceSymbol (const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var) const;
  void optimize (Node** parent) {}
  size_t stackCount(size_t sp, size_t& max_stack) { if (max_stack < sp+1) max_stack = sp+1; return sp; }
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
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const { out << "P["<<idx<<"]"; }
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent) {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const {}
  Node* derivative (const Node* var) const;
  void optimize (Node** parent) {}
  size_t stackCount(size_t sp, size_t& max_stack) { if (max_stack < sp+1) max_stack = sp+1; return sp; }
private:
  size_t idx;
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
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  // replaces the symbol with the node
  // parent points to the pointer of the current class
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  Node* derivative (const Node* var) const;
  void optimize (Node** parent);
  size_t stackCount(size_t sp, size_t& max_stack);
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
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  Node* derivative (const Node* var) const;
  // specific
  void addArgument (size_t pos, Node* arg, double m) { args[pos] = arg; mul[pos] = m; }
  void optimize (Node** parent);
  size_t stackCount(size_t sp, size_t& max_stack);
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
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  // needs to handle the divisions
  Node* derivative (const Node* var) const ;
  // specific
  void optimize (Node** parent);
  size_t stackCount(size_t sp, size_t& max_stack);
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
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  void print (std::ostream& out) const;
  bool operator != (const Node& node) const;
  bool operator < (const Node& node) const;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent);
  void replaceSymbol(const Node& node, Node** parent);
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const;
  Node* derivative (const Node* var) const;
  void optimize (Node** parent);
  size_t stackCount(size_t sp, size_t& max_stack);
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
  void derivative (Expression& ex, const Node* nd) const
  {
    if(root != nullptr)
    {
      Node* deri = root->derivative (nd);
      deri -> optimize (&deri);
      ex.fromNode (deri);
    }
  }
  void print (std::ostream& out) const { if (root != nullptr) root -> print (out); }
  Node* copy () const { if(root != nullptr) return root -> copy (); else return nullptr; }
  void stackCount (size_t& max_stack) const { if(root != nullptr) root -> stackCount (0, max_stack); }
  void evaluate (std::vector<Value>& stack, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
  {
    if (root != nullptr) root -> evaluate (stack, 0, var, par, len);
  }
  void evaluate (std::vector<double>& stack, const std::vector<double>& var, const std::vector<double>& par) const
  {
    if (root != nullptr) root -> evaluate (stack, 0, var, par);
  }
  void knutSplit (
    std::string& sysName,
    std::vector<std::string>& varName,
    std::vector<Expression>& varDotExpr,
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
