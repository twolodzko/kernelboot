#include <Rcpp.h>

class Foo
{
private:
  int x;
public:
  Foo() { x = 0; };
  void set(int y) { x = y; }
  int get() { return x; }
};


// [[Rcpp::export]]
void foo() {
  std::vector<Foo> xx(5);
  for (int i = 0; i < 5; i++)
    xx[i].set(i+1);
  for (int i = 0; i < 5; i++)
    Rcpp::Rcout << xx[i].get() << std::endl;
}

