
//: This Swift IOS playground shows how to solve a linear equation with n unknowns and n equations using guassian elimination.

//: Ax=b where A is a matrix with the coefficients of the linear equations, x and b are column vectors. x is the solution to the system of equations. The number of rows of x has to equal the number of columns of A from the rules of matrix multiplication.

//: ![System of equations](NumberedEquation2.gif)

//: Form the augmented matrix by adding a column to the end of A with the data from the column vector b.

//: ![augmented matrix for guassian elimination](augmentedMatrix.gif)

//: Put the augmented matrix in upper triangular form by performing row elimination.

//: ![upper triangular form](upperTriangularForm.gif)

//: First solve for xk and use back substitution to solve the remaining x variables.


import UIKit

var str = "Hello, playground"

//This struct stores a column vector.

public struct Vector {
    let rows :Int
    var grid=[Double]()
    init(rows: Int) {
        self.rows = rows
        grid = Array(count: rows, repeatedValue: 0.0)
    }
    func indexIsValid(row: Int) -> Bool {
        return row >= 0 && row < rows
    }
    subscript(row: Int) -> Double {
        get {
            assert(indexIsValid(row), "Index out of range")
            return grid[row]
        }
        set {
            assert(indexIsValid(row), "Index out of range")
            grid[row] = newValue
        }
    }
}

//This struct stores a matrix
public struct Matrix {
    let rows: Int, columns: Int
    var grid: [Double]
    
    init(rows: Int, columns: Int) {
        self.rows = rows
        self.columns = columns
        grid = Array(count: rows * columns, repeatedValue: 0.0)
    }
    
    func indexIsValidForRow(row: Int, column: Int) -> Bool {
        return row >= 0 && row < rows && column >= 0 && column < columns
    }
    
    subscript(row: Int, column: Int) -> Double {
        get {
            assert(indexIsValidForRow(row,column: column), "Index out of range")
            return grid[(row * columns) + column]
        }
        set {
            assert(indexIsValidForRow(row, column: column), "Index out of range")
            grid[(row * columns) + column] = newValue
        }
    }
}



//This is the gauss eliniation to solve a linear equation
func gauss(a:Matrix)->[Double]{
    //gauss elimination to solve linear equation. Use the augmented matrix!
    let rows=a.rows
    let columns=a.columns
    
    // matrix A will be the augmented matrix in upper triangular form
    var A=Matrix(rows: rows, columns: columns)
    A=a
    
    var n=rows
    var c: Double
    
    //x will be the solution in a Double array
    var x: [Double]
    x = Array(count: rows, repeatedValue: 0.0)
    
        // loop for the generation of upper triangular matrix
    /*
     
     This is the code in c
        
     for(j=1; j<=n; j++) /* loop for the generation of upper triangular matrix*/
     {
         for(i=1; i<=n; i++)
             {
             if(i>j)
             {
                 c=A[i][j]/A[j][j];
                 for(k=1; k<=n+1; k++)
                 {
                 A[i][k]=A[i][k]-c*A[j][k];
                 }
             }
         }
     }
     
     */
    // loop for the generation of upper triangular matrix
    //print("\(A)")
    for j in 0..<n{
        for i in 0..<n{
            if i>j{
            c=A[i,j]/A[j,j]
                for k in 0..<n+1{
                    A[i,k] -= c*A[j,k]
                    //uncomment this for trouble shooting
                    //print("A[\(i),\(k)]=\(A[i,k])")
                }
            }
        }
    }
    
    //print("\(A)")
    n -= 1
    x[n]=A[n,n+1]/A[n,n]
     //print("x[\(n)]=\(x[n])")

    var sum=0.0
     /* this loop is for backward substitution in c
     for(i=n-1; i>=1; i--)
     {
         sum=0;
             for(j=i+1; j<=n; j++)
             {
             sum=sum+A[i][j]*x[j];
             }
         x[i]=(A[i][n+1]-sum)/A[i][i];
     }
     */
    var i :Int

    for(i=n-1; i>=0; i -= 1){
        sum=0.0
        
        for j in i+1..<rows{
            sum += A[i,j]*x[j]
            //print("sum[\(j)]=\(sum)")
            
        }
        x[i]=(A[i,n+1]-sum)/A[i,i]
         //print("x[\(n)]=\(x[n])")
        
    }
    return x

}

func mvmul (a: Matrix, b: Vector) -> [Double]{
    //This function multiplies an M-by-P matrix A by a M-by-1 column vector B and stores the results in an M-by-1 vector C.
    
    let rows=a.rows
    let columns=a.columns
    let rows1=b.rows
    
    var a1=[Double]()
    var b1=[Double]()
    var count=0
    
    a1 = [Double](count:rows*columns, repeatedValue:0.0)
    for i in 0..<columns{
        for j in 0..<rows{
            print("\(count)")
            a1[count]=a[i,j]
            count += 1
        }
    }
    print("\(a1)")
    
    b1 = [Double](count:rows1, repeatedValue:0.0)
    
    count=0
    for i in 0..<rows1{
        print("\(count)")
        b1[count]=b[i]
        count += 1
    }
    print("\(b1)")
    
    var result = [Double](count:rows1, repeatedValue:0.0)
    
    for i in 0..<columns{
        for j in 0..<rows1{
            result[i] += a[i,j]*b[j]
        }
    }
    
    return result
}
/********************************
Now check the solution by using matrix multiplication to check the solution x
 ************************************/

var m1=Matrix(rows: 3,columns: 4)

m1[0, 0] = 10.0
m1[0, 1] = -7.0
m1[0, 2] = 3.0
m1[0, 3] = 5.0
m1[1, 0] = -6.0
m1[1, 1] = 8.0
m1[1, 2] = 4.0
m1[1, 3] = 7.0
m1[2, 0] = 2.0
m1[2, 1] = 6.0
m1[2, 2] = 9.0
m1[2, 3] = -1.0

print("\(m1)")


let solution=gauss(m1)

print("\(solution)")

//check

/********************************
 Now check the solution by using matrix multiplication to check the solution x
 ************************************/

var solution1=Vector(rows: 3)

solution1[0] = -7.809086
solution1[1] = -8.690904
solution1[2] = 7.418178


var m2=Matrix(rows: 3,columns: 3)

m2[0, 0] = 10.0
m2[0, 1] = -7.0
m2[0, 2] = 3.0
m2[1, 0] = -6.0
m2[1, 1] = 8.0
m2[1, 2] = 4.0
m2[2, 0] = 2.0
m2[2, 1] = 6.0
m2[2, 2] = 9.0



let checkB=mvmul(m2, b: solution1)

print("\(checkB)")

