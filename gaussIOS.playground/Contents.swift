


//: This Swift IOS playground shows how to solve a linear equation with n unknowns and n equations using guassian elimination.

//: Ax=b where A is a matrix with the coefficients of the linear equations, x and b are column vectors. x is the solution to the system of equations. The number of rows of x has to equal the number of columns of A from the rules of matrix multiplication.

//: ![System of equations](NumberedEquation2.gif)

//: Form the augmented matrix by adding a column to the end of A with the data from the column vector b.

//: ![augmented matrix for guassian elimination](augmentedMatrix.gif)

//: Put the augmented matrix in upper triangular form by performing row elimination.

//: ![upper triangular form](upperTriangularForm.gif)

//: First solve for xk and use back substitution to solve the remaining x variables.


import UIKit


//left this statement in to verify playground is running
var str = "Hello, playground"

//: ### struct column vector

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

//: ### struct for matrix

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

//: ### gauss function for solving n equations in n unknowns

//This is the gauss eliniation to solve a linear equation
func gauss(a:Matrix)->Vector{
    //gauss elimination to solve linear equation. Use the augmented matrix!
    let rows=a.rows
    let columns=a.columns
    
    // matrix A will be the augmented matrix in upper triangular form
    var A=Matrix(rows: rows, columns: columns)
    A=a
    
    var n=rows
    var c: Double
    
    //x will be the solution in a Double array
    
    var x=Vector(rows: rows)
    
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
    var i=n
    //print("n=\(n)")
    //for(i=n-1; i>=0; i -= 1){
    //eliminated c type for statement per Swift warning
    for _ in 0..<n{
        i -= 1

        //print("i=\(i)")
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

//: ### function mvmul matrix times vector multiplication

func mvmul (a: Matrix, b: Vector) -> Vector{
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
            //print("\(count)")
            a1[count]=a[i,j]
            count += 1
        }
    }
    print("\(a1)")
    
    b1 = [Double](count:rows1, repeatedValue:0.0)
    
    count=0
    for i in 0..<rows1{
        //print("\(count)")
        b1[count]=b[i]
        count += 1
    }
    print("\(b1)")
    
    var result=Vector(rows: rows)
    
    for i in 0..<columns{
        for j in 0..<rows1{
            result[i] += a[i,j]*b[j]
        }
    }
    
    return result
}
/********************************
Lets try to solve the first system of equations Ax=b
 ************************************/

//: Lets try to solve the first system of equations Ax=b

//: ![first system of equtions to solve](guass1.gif)

//: First put the system of equatons in augmented form

//: ![first system of equtions to solve](guass2.gif)
var m0=Matrix(rows: 3,columns: 4)

m0[0, 0] = 9.0
m0[0, 1] = 3.0
m0[0, 2] = 4.0
m0[0, 3] = 7.0
m0[1, 0] = 4.0
m0[1, 1] = 3.0
m0[1, 2] = 4.0
m0[1, 3] = 8.0
m0[2, 0] = 1.0
m0[2, 1] = 1.0
m0[2, 2] = 1.0
m0[2, 3] = 3.0

print("Augmented matrix \(m0)")


let solution=gauss(m0)

print("The solution to the first example is \(solution)")

//check

/********************************
 Now check the solution by using matrix multiplication to check the solution x
 ie. Make sure matrix A times vector x solution equals b
 ************************************/


//: Now check the solution by using matrix multiplication to check the solution x
//: ie. Make sure matrix A times solution vector x equals column vector b


var m1=Matrix(rows: 3,columns: 3)

m1[0, 0] = 9.0
m1[0, 1] = 3.0
m1[0, 2] = 4.0
m1[1, 0] = 4.0
m1[1, 1] = 3.0
m1[1, 2] = 4.0
m1[2, 0] = 1.0
m1[2, 1] = 1.0
m1[2, 2] = 1.0

print("matrix A \(m1)")


let checkB=mvmul(m1, b: solution)

print("This is the check for b \(checkB)")

//: ### solution checks for b
var b=Vector(rows: 3)
b[0] = 7.0
b[1] = 8.0
b[2] = 3.0

func areVectorsEqual (a: Vector, b: Vector) -> Bool {
    //This function returns true is 2 input vectors are equal
    
    assert(a.rows == b.rows, "Expected vectors of the same length, instead got arrays of two different lengths")
    
    var result = false
    for index in 0..<a.rows {
        if a[index]==b[index]{result=true}
        else{result=false
            return result
        }
    }
    return result
}

if areVectorsEqual(b, b: checkB){
    print("b and checkB vectors are equal\n gaussian elimination functions checks Ax=b")
}
else{print("solution does not check/n b and checkB vectors are not equal")}




