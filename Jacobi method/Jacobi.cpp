#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

class ColumnVector {

private:

    int size=0;
    double *v=NULL;

public:

    ColumnVector() {}

    ColumnVector(int size) {
        this->size=size;
        const int param = size;
        this->v = new double[param];
    }

    ColumnVector(int size,double a[]) {
        this->size=size;
        const int param = size;
        v = new double[param];
        for (int i=0;i<size;i++) {
            this->v[i]=a[i];
        }
    }

    int getSize() const {
        return size;
    }

    double get(int index) const {
        return this->v[index];
    }

    void set(int index,double val) {
        this->v[index]=val;
    }

    double get_length() const {
        double sum=0;
        for (int i=0;i<size;i++) {
            sum+=v[i]*v[i];
        }
        sum=sqrt(sum);
        return sum;
    }

    ColumnVector& operator=(const ColumnVector &V) {
        this->size=V.getSize();
        v = new double[size];
        for (int i=0;i<size;i++) {
            this->v[i]=V.get(i);
        }
        return *this;
    }

    ColumnVector operator+(const ColumnVector &V) {
        ColumnVector sum(V.getSize());
        for (int i=0;i<size;i++) {
            sum.set(i,this->v[i]+V.get(i));
        }
        return sum;
    }

    ColumnVector operator-(const ColumnVector &V) {
        ColumnVector diff(V.getSize());
        for (int i=0;i<size;i++) {
            diff.set(i,this->v[i]-V.get(i));
        }
        return diff;
    }

    friend ostream & operator << (ostream &out, const ColumnVector &V){
        for (int i=0;i<V.size;i++) {
            out<<V.v[i]<<'\n';
        }
        return out;
    }

    friend istream & operator >> (istream &in, ColumnVector &V){
        in>>V.size;
        V.v = new double[V.size];
        for (int i=0;i<V.size;i++) {
            in>>V.v[i];
        }
        return in;
    }
};

class Matrix {
protected:
    int rows=0,columns=0;
    double **m=NULL;
public:
    Matrix() {}

    Matrix(int rows,int columns) {
        this->rows=rows;
        this->columns=columns;
        m = new double*[rows];
        for (int i=0;i<rows;i++) {
            m[i] = new double[columns];
        }
        for (int i=0;i<rows;i++) {
            for (int j=0;j<columns;j++) {
                m[i][j]=0;
            }
        }
    }

    Matrix(Matrix A,Matrix B) {
        this->rows=A.getRows();
        this->columns=A.getColumns()+B.getColumns();
        m = new double*[rows];
        for (int i=0;i<rows;i++) {
            m[i] = new double[columns];
        }
        for (int i=0;i<rows;i++) {
            for (int j=0;j<columns;j++) {
                if (j<A.getColumns()) m[i][j]=A.get(i,j);
                else m[i][j]=B.get(i,j-A.getColumns());
            }
        }
    }

    int getRows() const {
        return rows;
    }

    int getColumns() const {
        return columns;
    }

    void set(int i,int j,double val) {
        m[i][j]=val;
    }

    double get(int i,int j) const{
        return m[i][j];
    }

    Matrix operator+(const Matrix &M) {
        if (M.getRows()!=rows || M.getColumns()!=columns) {
            cout<<"Error: the dimensional problem occurred\n";
            return *(new Matrix());
        }
        else {
            Matrix sum(rows,columns);
            for (int i=0;i<rows;i++) {
                for (int j=0;j<columns;j++) {
                    sum.set(i,j,this->m[i][j]+M.get(i,j));
                }
            }
            return sum;
        }
    }

    Matrix operator-(const Matrix &M) {
        if (M.getRows()!=rows || M.getColumns()!=columns) {
            cout<<"Error: the dimensional problem occurred\n";
            return *(new Matrix());
        }
        else {
            Matrix dif(rows,columns);
            for (int i=0;i<rows;i++) {
                for (int j=0;j<columns;j++) {
                    dif.set(i,j,this->m[i][j]-M.get(i,j));
                }
            }
            return dif;
        }
    }

    Matrix operator*(const Matrix &M) {
        if (columns!=M.getRows()) {
            cout<<"Error: the dimensional problem occurred\n";
            return *(new Matrix());
        }
        else {
            Matrix product(rows,M.getColumns());
            for (int i=0;i<rows;i++) {
                for (int j=0;j<M.getColumns();j++) {
                    double sum=0;
                    for (int k=0;k<columns;k++) {
                        sum+=m[i][k]*M.get(k,j);
                    }
                    product.set(i,j,sum);
                }
            }
            return product;
        }
    }

    ColumnVector operator*(const ColumnVector &V) {
        if (columns!=V.getSize()) {
            cout<<"Error: the dimensional problem occurred\n";
            return *(new ColumnVector());
        }
        else {
            ColumnVector product(rows);
            for (int i=0;i<rows;i++) {
                double sum=0;
                for (int k=0;k<columns;k++) {
                    sum+=m[i][k]*V.get(k);
                }
                product.set(i,sum);
            }
            return product;
        }
    }

    Matrix& operator=(const Matrix &M) {
        this->rows=M.getRows();
        this->columns=M.getColumns();
        m = new double*[rows];
        for (int i=0;i<rows;i++) {
            m[i] = new double[columns];
        }
        for (int i=0;i<rows;i++) {
            for (int j=0;j<columns;j++) {
                m[i][j]=M.get(i,j);
            }
        }
        return *this;
    }

    void transpose() {
        double **temp = new double*[columns];
        for (int j=0;j<columns;j++) temp[j]=new double[rows];
        for (int i=0;i<rows;i++) {
            for (int j=0;j<columns;j++) {
                temp[j][i]=m[i][j];
            }
        }
        m = new double*[columns];
        for (int j=0;j<columns;j++) m[j]=new double[rows];
        swap(rows,columns);
        for (int i=0;i<rows;i++) {
            for (int j=0;j<columns;j++) {
                m[i][j]=temp[i][j];
            }
        }
    }

    friend ostream & operator << (ostream &out, const Matrix &M){
        for (int i=0;i<M.rows;i++) {
            for (int j=0;j<M.columns;j++) {
                out<<M.m[i][j]<<' ';
            }
            out<<endl;
        }
        return out;
    }

    friend istream & operator >> (istream &in, Matrix &M){
        in>>M.rows>>M.columns;
        M.m = new double*[M.rows];
        for (int i=0;i<M.rows;i++) {
            M.m[i] = new double[M.columns];
        }
        for (int i=0;i<M.rows;i++) {
            for (int j=0;j<M.columns;j++) {
                in>>M.m[i][j];
            }
        }
        return in;
    }
};

class SquareMatrix : public Matrix {

protected:
    int n;

public:
    SquareMatrix() : Matrix() {}

    SquareMatrix(int n) : Matrix(n,n) {this->n=n;}

    SquareMatrix(Matrix M) : Matrix (M.getRows(),M.getColumns()) {
        this->n=M.getRows();
        this->rows=n;
        this->columns=n;
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                this->m[i][j]=M.get(i,j);
            }
        }
    }

    int getN() const{
        return n;
    }

    SquareMatrix& operator=(const Matrix &M) {
        this->rows=M.getRows();
        this->columns=M.getColumns();
        this->n=M.getRows();
        m = new double*[rows];
        for (int i=0;i<rows;i++) {
            m[i] = new double[columns];
        }
        for (int i=0;i<rows;i++) {
            for (int j=0;j<columns;j++) {
                m[i][j]=M.get(i,j);
            }
        }
        return *this;
    }

    SquareMatrix operator-(const SquareMatrix &M) {
        if (M.getRows()!=rows || M.getColumns()!=columns) {
            cout<<"Error: the dimensional problem occurred\n";
            return *(new SquareMatrix());
        }
        else {
            SquareMatrix dif(rows);
            for (int i=0;i<rows;i++) {
                for (int j=0;j<columns;j++) {
                    dif.set(i,j,this->m[i][j]-M.get(i,j));
                }
            }
            return dif;
        }
    }

    friend ostream & operator << (ostream &out, const SquareMatrix &M){
        for (int i=0;i<M.n;i++) {
            for (int j=0;j<M.n;j++) {
                out<<M.m[i][j]<<' ';
            }
            out<<endl;
        }
        return out;
    }

    friend istream & operator >> (istream &in, SquareMatrix &M){
        in>>M.n;
        M.rows=M.n;
        M.columns=M.n;
        M.m = new double*[M.n];
        for (int i=0;i<M.n;i++) {
            M.m[i] = new double[M.n];
        }
        for (int i=0;i<M.n;i++) {
            for (int j=0;j<M.n;j++) {
                in>>M.m[i][j];
            }
        }
        return in;
    }
};

class IdentityMatrix : public SquareMatrix {

public:
    IdentityMatrix() : SquareMatrix() {}

    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i=0;i<n;i++) m[i][i]=1;
    }

};

class EliminationMatrix : public IdentityMatrix {
private:
    int i,j;
    Matrix x;

public:
    EliminationMatrix() : IdentityMatrix() {}

    EliminationMatrix(int i,int j,Matrix x) : IdentityMatrix(x.getRows()) {
        this->i=i;
        this->j=j;
        this->x=x;
        m[i][j]=-x.get(i,j)/x.get(j,j);
    }
};

class PermutationMatrix : public IdentityMatrix {

private:
    int i,j;

public:

    PermutationMatrix() : IdentityMatrix() {}

    PermutationMatrix(int i,int j,int n) : IdentityMatrix(n) {
        this->i=i;
        this->j=j;
        swap(m[i][j],m[i][i]);
        swap(m[j][i],m[j][j]);
    }
};

SquareMatrix inv(SquareMatrix A) {
    int n=A.getN();
    IdentityMatrix I(n);
    SquareMatrix inverse = I;
    int x=1,y=0;
    int step=1;
    for (int i=0;i<((n)*(n-1))/2;i++) {
        double piv=A.get(y,y);
        int s=-1;
        for (int j=x;j<n;j++) {
            if (abs(A.get(j,y))>abs(piv)) {
                piv=A.get(j,y);
                s=j;
            }
        }
        if (s!=-1) {
            PermutationMatrix P(y,s,n);
            A=P*A;
            inverse=P*inverse;
            step++;
        }
        EliminationMatrix E(x,y,A);
        bool skip=(A.get(x,y)==0);
        x++;
        if (x==n) {y++; x=y+1;}
        if (skip) continue;
        step++;
        A=E*A;
        inverse=E*inverse;
    }
    y=n-1;
    x=y-1;
    for (int i=0;i<((n)*(n-1))/2;i++) {
        EliminationMatrix E(x,y,A);
        bool skip=(A.get(x,y)==0);
        x--;
        if (x==-1) {y--; x=y-1;}
        if (skip) continue;
        step++;
        A=E*A;
        inverse=E*inverse;
    }
    for (int i=0;i<n;i++) {
        double d=A.get(i,i);
        for (int j=0;j<A.getColumns();j++) {
            if (abs(A.get(i,j))>0.000001 && d) A.set(i,j,A.get(i,j)/d);
            if (abs(inverse.get(i,j))>0.000001 && d) inverse.set(i,j,inverse.get(i,j)/d);
        }
    }
    return inverse;
}

double det(SquareMatrix A) {
    int n=A.getN();
    int x=1,y=0;
    int step=1;
    for (int i=0;i<((n)*(n-1))/2;i++) {
        int piv=A.get(y,y);
        int s=-1;
        for (int j=x;j<n;j++) {
            if (abs(A.get(j,y))>abs(piv)) {
                piv=A.get(j,y);
                s=j;
            }
        }
        if (s!=-1) {
            PermutationMatrix P(y,s,n);
            A=P*A;
            step++;
        }
        EliminationMatrix E(x,y,A);
        bool skip=(A.get(x,y)==0);
        x++;
        if (x==n) {y++; x=y+1;}
        if (skip) continue;
        step++;
        A=E*A;
    }
    double res=1;
    for (int i=0;i<n;i++) {
        res*=A.get(i,i);
    }
    return res;
}

bool notDomin(SquareMatrix A) {
    for (int i=0;i<A.getN();i++) {
        double sum=0;
        for (int j=0;j<A.getN();j++) {
            if (i==j) continue;
            sum+=abs(A.get(i,j));
        }
        if (abs(A.get(i,i))<sum) return true;
    }
    return false;
}

double power(double x,int n) {
    double p=1;
    for (int i=0;i<n;i++) p*=x;
    return p;
}

double x[100],y[100];

int main()
{
    cout<<fixed<<setprecision(4);
    int m;
    cin>>m;
    for (int i=0;i<m;i++) {
        cin>>x[i]>>y[i];
    }
    int n;
    cin>>n;
    Matrix A(m,n+1);
    for (int i=0;i<m;i++) {
        for (int j=0;j<n+1;j++) {
            A.set(i,j,power(x[i],j));
        }
    }
    ColumnVector b(m,y);
    A.transpose();
    Matrix A_T=A;
    A.transpose();
    cout<<"A:\n"<<A;
    cout<<"A_T*A:\n"<<A_T*A;
    cout<<"(A_T*A)^-1:\n"<<inv(SquareMatrix(A_T*A));
    cout<<"A_T*b:\n"<<A_T*b;
    cout<<"x~:\n";
    ColumnVector x= inv(A_T*A)*A_T*b;
    cout<<x;
}


