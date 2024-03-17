//Newton Raphson Helper Functions

function createPoly(arr){//turns nx2 array into degree n complex polynomial

    let P = [];
    for(let coeff of arr){
        P.push(new Complex(coeff[0], coeff[1]));
    }

    return P;
}

function evaluate(P, z){ //evaluates polynomial P at z

    let s = P[0];
    let z0 = new Complex(1,0);

    for(let i=1; i<P.length; i++){
        z0 = z0.mult(z);
        s=s.add(P[i].mult(z0));
    }

    return s;
}

function derivative(P){ //returns derivative of polynomial

    if(P.length==1){
        return [new Complex(0,0)];
    }

    let dP = [];
    for(let i=1; i<P.length; i++){
        dP.push(P[i].mult(new Complex(i, 0)));
    }

    return dP;

}

function newton_raphson(P, z0, niter){
    let z = z0;
    let dP = derivative(P);
    let q = rational(P, dP);

    for(let i=0; i<niter; i++){
        z = z.sub(q(z));
    }

    return z;
}

function remove_trailing_zeros(P){

    let n = 0;
    for(let i=P.length-1; i>0; i--){
        if(P[i].abs > 0.000000001){
            n = i;
            break;
        }
    }

    return jsn.parse(P, 0, n+1, 1);
}

function same_length(P,Q){
    
    let n = P.length;
    let m = Q.length;

    if(n<m){
        for(let i=0; i<m-n; i++){
            P.push(new Complex(0,0));
        }
    }
    else if(n>m){
        for(let i=0; i<n-m; i++){
            Q.push(new Complex(0,0));
        }
    }
}

function multiply_poly(P,Q){ //multiply P and Q

    same_length(P, Q);

    let vP = jsn.fft(jsn.fpad(P, 1));
    let vQ = jsn.fft(jsn.fpad(Q, 1));

    let prod = [];
    for(let i=0; i<vP.length; i++){
        prod.push(vP[i].mult(vQ[i]));
    }

    return remove_trailing_zeros(jsn.ifft(prod));
}

function subtract_poly(P, Q){

    same_length(P, Q);
    let h = [];

    for(let i=0; i<P.length; i++){
        h.push(P[i].sub(Q[i]));
    }

    return h;
}

function divide_poly(P, Q){ //assuming Q is a degree 1 polynomial which divides P

    same_length(P, Q);

    let vP = jsn.fft(jsn.fpad(P));
    let vQ = jsn.fft(jsn.fpad(Q));

    //locate zero of Q
    let n = -1;
    for(let i=0; i<vQ.length; i++){
        if(vQ[i].abs < 0.0000000000000001){
            n = i;
            break;
        }
    }

    //divide values

    let h = [];

    let m = vP.length;
    let woff = Complex.polar(1, 2*Math.PI*n / m).add(new Complex(0.00001, 0.00001));

    for(let i=0; i<m; i++){
        if(i==n){
            h.push(evaluate(P, woff).div(evaluate(Q, woff)));
        }
        else{
            h.push(vP[i].div(vQ[i]));
        }
    }

    return remove_trailing_zeros(jsn.ifft(h));

}

function rational(P, Q){ //creates rational function from polynomials, represented as coefficient arrays
    
    return (z)=>{

        if(z.re == Infinity){
            if(P.length>Q.length){ //this line here assumes that the polynomials are represented without trailing zeros
                return z;
            }
            else if(P.length<Q.length){
                return new Complex(0,0);
            }
            else{
                return P[P.length-1].div(Q[Q.length-1]);
            }
        }

        else{

            let pz = evaluate(P, z);
            let qz = evaluate(Q, z);

            if(qz.abs<0.000000001){
                return new Complex(Infinity, Infinity);
            }

            else{
                return pz.div(qz);
            }
        }
    }

}


function findRoots(P){

    let n = P.length-1;

    let roots = [];

    for(let i=0; i<n; i++){

        //generate guess for newton raphson method
        let x = jsn.random(-1, 1);
        let y = jsn.random(-1, 1);
    
        let z0 = new Complex(x,y);

        let r = newton_raphson(P, z0, 20);
        roots.push(r);

        P = divide_poly(P, [r.mult(new Complex(-1,0)), new Complex(1,0)]);
    }

    return roots;
}


function genady_polynomial(c1, p, q){

    let genady = [new Complex(1,0)];

    //start by calculating the inner summands, and keep them in an array
    let derv = new Complex(1,0).div(c1.mult(new Complex(2,0)))
    let z0 = c1;


    let derivs = [derv];

    for(let i=2; i<p*q+1; i++){
        z0 = z0.mult(z0).add(c1);
        if(derv.abs<0.00000000000001){
            derv = new Complex(0,0);
        }
        else{
            derv = derv.div(new Complex(2,0)).div(z0);
        }
        derivs.push(derv);
    }
    //now calculate polynomial
    for(let n=1; n<q; n+=1){
        genady.push(jsn.parse(derivs, p*(n-1)+1, p*n+1, 1).reduce((sum, a)=> sum.add(a), new Complex(0,0)));
    }

    genady.push(jsn.parse(derivs, p*(q-1)+1, p*q, 1).reduce((sum, a)=> sum.add(a), new Complex(0,0)));

    return genady;

}


//calculate all roots 

let c1 = new Complex(-1.4, -0.01);
let p = 3;
let q = 2*2*2*2*2;

let g = genady_polynomial(c1, p , q);
console.log(findRoots(g));