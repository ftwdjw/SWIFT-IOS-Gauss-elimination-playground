


//: This Swift IOS playground shows how to solve a linear equation with n unknowns and n equations using guassian elimination.

//: Ax=b where A is a matrix with the coefficients of the linear equations, x and b are column vectors. x is the solution to the system of equations. The number of rows of x has to equal the number of columns of A from the rules of matrix multiplication.

//: ![System of equations](NumberedEquation2.gif)

//: Form the augmented matrix by adding a column to the end of A with the data from the column vector b.

//: ![augmented matrix for guassian elimination](augmentedMatrix.gif)

//: Put the augmented matrix in upper triangular form by performing row elimination.

//: ![upper triangular form](upperTriangularForm.gif)

//: First solve for xk and use back substitution to solve the remaining x variables.


import UIKit

//: ### This struct stores a column vector.

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

//: ### This struct stores a matrix.

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

//: ### gauss function for solving a system of equations n equations in n unknowns Ax=b

//This is the gauss eliniation to solve a linear equation
func gauss(a:Matrix, b:Vector)->(x: Vector,isValid: String){
        assert(a.rows == b.rows, "The number of rows of the matrix has to be equal to the number of rows of the b vector")
    //gauss elimination to solve linear equation Ax=b. inputs:a is the A matrix. b is the b vector. outputs: x is the solution vector. isValid is a string to indicate if the result is valid. Some systems of equations have no solution.
    let rows=a.rows
    let columns=a.columns
    
    // matrix A will be the augmented matrix in upper triangular form
    var A=Matrix(rows: rows, columns: columns+1)
    //print("\(A)")
    
    //fill in augmented matrix with matrix A and vector b in the last column
    for i in 0..<rows{
        for j in 0..<columns+1{
            if j<columns{
                 A[i,j]=a[i,j]
            }
            else{
                A[i,j]=b[i]
            }
      
        }}
     //print("\(A)")
    
    var n=rows
    var c: Double
    
    //x will be the solution in a Vector
    
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
    
    //print("This is the upper triangular matrix\(A)")
    
    n -= 1
    //added assertion to make sure there is a solution for the system of equations
    //assert(A[n,n] != 0.0, "No solution for this system of equations!")
    
    if A[n,n]==0 {
        //print out x=0 and no valid solution if there is no solution
        return (x: x, isValid: "There is a no valid solution")
    }
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
    return (x: x, isValid: "There is a valid solution")

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
    //print("\(a1)")
    
    b1 = [Double](count:rows1, repeatedValue:0.0)
    
    count=0
    for i in 0..<rows1{
        //print("\(count)")
        b1[count]=b[i]
        count += 1
    }
    //print("\(b1)")
    
    var result=Vector(rows: rows)
    
    for i in 0..<columns{
        for j in 0..<rows1{
            result[i] += a[i,j]*b[j]
        }
    }
    
    return result
}

//: ### function printEquation to print Ax=b is row and column form

func printEquation (A:Matrix, b:Vector){
    let columns=A.columns
    let rows=A.rows
    //var rowToPrint=[Double]()
    //rowToPrint.append(Array(count: rows, repeatedValue:0.0))
    print("\n")
    var count=0
    print("Ax=b in row form (matrix equation)\n")
    for i in 0..<columns{
         print("[", terminator:"")
        for j in 0..<rows{
            print("A[\(i+1),\(j+1)])=\(A[i,j]) ", terminator:"")
        }
        print("]", terminator:"")
         print("[", terminator:"")
        print(" x[\(count+1)] ", terminator:"")
        print("]", terminator:"")
        count += 1
        print(" = [b[\(i+1)]=\(b[i]) ", terminator:"")
        print("]", terminator:"")
        print("\n")
    }
    
    count=0
    
    print("Ax=b in column form (vector equation)\n")
    for i in 0..<columns{
        for j in 0..<rows{
            if i==0{
                if count==0{
                    print("x[\(count+1)] [A[\(i+1),\(j+1)])=\(A[i,j])]  ", terminator:"")}
                else{
                    print("+x[\(count+1)][A[\(i+1),\(j+1)])=\(A[i,j])]  ", terminator:"")
                }
                
                
                
                }
            else{
                print("     [A[\(i+1),\(j+1)])=\(A[i,j])]  ", terminator:"")
            }
            
            count += 1
            if count==columns{print("= [b[\(i+1)]=\(b[i])] ", terminator:"")}
        }
        print("\n")
        count=0
    }
    
}//end print func

//: ### function printSolution

func printSolution(a:Vector){
    let rows=a.rows
    var solution=Vector(rows: rows)
    solution=a
    
    print("The solution to the system of equations is:")
    
    for index in 0..<rows{
        print("x[\(index+1)]=\(solution[index])")
    }
}

/********************************
Lets try to solve the first system of equations Ax=b
 ************************************/

//: Lets try to solve the first system of equations Ax=b

//: ![first system of equtions to solve](guass1.gif)

//: First put the system of equatons in augmented form

//: ![first system of equtions to solve](guass2.gif)

var m0=Matrix(rows: 3,columns: 3)

m0[0, 0] = 9.0
m0[0, 1] = 3.0
m0[0, 2] = 4.0
m0[1, 0] = 4.0
m0[1, 1] = 3.0
m0[1, 2] = 4.0
m0[2, 0] = 1.0
m0[2, 1] = 1.0
m0[2, 2] = 1.0

print("Matrix A \(m0)")

var b=Vector(rows: 3)
b[0] = 7.0
b[1] = 8.0
b[2] = 3.0

print("Vector b \(b)")

printEquation(m0,b:b)


let solution=gauss(m0, b: b).x

printSolution(solution)

//check

/********************************
 Now check the solution by using matrix multiplication to check the solution x
 ie. Make sure matrix A times vector x solution equals b
 ************************************/


//: Now check the solution by using matrix multiplication to check the solution x
//: ie. Make sure matrix A times solution vector x equals column vector b




let checkB=mvmul(m0, b: solution)

print("\n\nThis is the check for b \(checkB)")



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

//: Now show an example with no solution
print("\n\nNow show an example with no solution")

var m2=Matrix(rows: 2,columns: 2)

m2[0, 0] = 1.0
m2[0, 1] = -1.0
m2[1, 0] = 2.0
m2[1, 1] = -2.0

print("Matrix A \(m2)")

var b1=Vector(rows: 2)
b1[0] = 4.0
b1[1] = -4.0

print("b Vector \(b1)")

printEquation(m2,b:b1)

let NoSolution=gauss(m2,b: b1).x
let NoSolutionString=gauss(m2,b: b1).isValid

print( NoSolutionString)





