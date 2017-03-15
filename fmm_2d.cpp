#include<iostream>
#include<complex>
#include<cmath>
#include<queue>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
#define MAX_SIZE 341//size of Q which contains the nodes at same level;currently 6 levels are taken

class node{
   public:
      float Q;
      int noc;
      int charge_index[10];  //assuming not more than 10 charges will be in a box at finest level
      int level;
      complex<double> loc;
      complex<double> Zo;
      complex<double> ak[100];
      complex<double> bk[100];//multipole expansion coeffecients  ;p < 100
      complex<double> local_ex_bk[100];//local expansion coeffecients ;;-ir
      
      queue<node*> nn;
      queue<node*> il;

      node *first_child;
      node *second_child;
      node *third_child;
      node *fourth_child;
      node *parent;
      
      /*//constructor
      node(){ 
      first_child = NULL;
      second_child = NULL;
      third_child = NULL;
      fourth_child = NULL;
      parent = NULL;
      noc = 0;
      Q = 0;
      }*/
};


class Queue_CA
{
public:
	node* A[MAX_SIZE];
	int front, rear; 
public:
	// Constructor - set front and rear as -1. 
	// We are assuming that for an empty Queue, both front and rear will be -1.
	Queue_CA()
	{
		front = -1; 
		rear = -1;
	}

	// To check wheter Queue is empty or not
	bool IsEmpty()
	{
		return (front == -1 && rear == -1); 
	}

	// To check whether Queue is full or not
	bool IsFull()
	{
		return (rear+1)%MAX_SIZE == front ? true : false;
	}
       
        int Size()
        {
           if(IsEmpty()) 
              return(0);
           else
              return((rear+MAX_SIZE-front)%MAX_SIZE + 1);
        }

	// Inserts an element in queue at rear end
	void Enqueue(node* x)
	{
		cout<<"Enqueuing "<<x<<" \n";
		if(IsFull())
		{
			cout<<"Error: Queue is Full\n";
			return;
		}
		if (IsEmpty())
		{ 
			front = rear = 0; 
		}
		else
		{
		    rear = (rear+1)%MAX_SIZE;
		}
		A[rear] = x;
	}

	// Removes an element in Queue from front end. 
	void Dequeue()
	{
		cout<<"Dequeuing \n";
		if(IsEmpty())
		{
			cout<<"Error: Queue is Empty\n";
			return;
		}
		else if(front == rear ) 
		{
			rear = front = -1;
		}
		else
		{
			front = (front+1)%MAX_SIZE;
		}
	}
	// Returns element at front of queue. 
	node* Front()
	{
		if(front == -1)
		{
			cout<<"Error: cannot return front from empty queue\n";
			return NULL; 
		}
		return A[front];
	}
	/* 
	   Printing the elements in queue from front to rear. 
	   This function is only to test the code. 
	   This is not a standard function for Queue implementation. 
	*/
	void Print()
	{
		// Finding number of elements in queue  
		int count = (rear+MAX_SIZE-front)%MAX_SIZE + 1;
		cout<<"Queue       : ";
		for(int i = 0; i <count; i++)
		{
			int index = (front+i) % MAX_SIZE; // Index of element while travesing circularly from front
			cout<<A[index]<<" ";
		}
		cout<<"\n\n";
	}
};


int preorder(node* node1,int l,double x,double y,complex<double>);
node* newnode(int l);
void levelorder(node* root);
void step3func(Queue_CA Q,int l);
void step3_il(Queue_CA Q,int l);
void step3_lexp(Queue_CA Q,int l);
complex<double> binom_coeff(int n,int r);
int fact(int r);
int sum_gp(int b);

int n,noc;
double p;
double a=16.000; //global variable; a= length of side of square at level 0
float q[100]; //assuming max 100 charges are there in the computational space
VectorXd qx,qy;
complex<double> qz[100];
complex<double> q_complex[100];


int main(){
    int j,k;
    double epsilon = 10^-16;
    int c = 2;
    p = (log(1/epsilon))/log(c);
    node *root = NULL;
    cout << "enter the number of levels : ";
    cin >> n;
    cout << "\nthe computational space is taken in the xy plane with x=[0,16) and y=[0,16)";
    cout << "\nenter the number of charges : ";
    cin >> noc;
    VectorXd q=5*(VectorXd::Random(noc)); //range: [-5,5]
    qx=16*((VectorXd::Ones(noc))+(VectorXd::Random(noc)))/2; //range: [0,16] location of charges in xy plane
    qy=16*((VectorXd::Ones(noc))+(VectorXd::Random(noc)))/2; 
    for(j=0;j<noc;j++){
        real(qz[j])=qx[j];
        imag(qz[j])=qy[j];
        real(q_complex[j]) = q[j];
    }
    preorder(root,0,a/2.0,a/2.0,0+0I);//forming tree and evaluating multipole expansion coeffecients
    levelorder(root);//evaluating local expansions
}

node* newnode(int l){
      int i;
      node* node1 = (node*)malloc(sizeof(node));
      node1->level = l;
      node1->first_child = NULL;
      node1->second_child = NULL;
      node1->third_child = NULL;
      node1->fourth_child = NULL;
      node1->noc = 0;
      node1->Q = 0;
      for(i=0;i<p;i++)
         node1->local_ex_bk[i] = 0+0I;
      return(node1);
}

//step 1 & 2
//creates tree datastrucrture for the problem;calculates the location of center for each box at the finest level;finds which charges are there in which box --done only at the finest level boxes

int preorder(node* node1,int l,double x,double y,complex<double> parent_loc){
    complex<double> complex_f;
    complex<double> sec_term1,sec_term2,sec_term3,sec_term4;
    static int j;//initialised to zero by default
    cout<<l<<"\n";
    j++;
    int d,i,k,r,f;
    if(node1 == NULL && l >= n+1){
      return(n);//returns parent level
     }
     else if(node1 == NULL){
          node1 = newnode(l);//now node1 is not null
          real(node1->loc) = x;
          imag(node1->loc) = y;         
          if(l == n){
            for(i=0;i<noc;i++){
               if((real(node1->loc)-a/(2^(l+1) < qx[i]) && (real(node1->loc)+a/(2^(l+1)) > qx[i]) && (imag(node1->loc)-a/(2^(l+1)) < qy[i]) && (imag(node1->loc))+a/(2^(l+1)) < qy[i])){
                  node1->noc = node1->noc + 1;
                  node1->charge_index[noc-1] = i;
                  node1->Q = node1->Q + q[i]; //multipole Q
                  node1->Zo = node1->loc-parent_loc;
               }
            }
          
            //step 1
            //multipole ak calculation starts from index 1
            for(f=1;f<=p;f++){
               real(complex_f) = f;
               for(d=1;d<=node1->noc;d++){
                    node1->ak[f] = node1->ak[f] + q_complex[node1->charge_index[d]]*pow(qz[node1->charge_index[d]]-(node1->loc),f)/complex_f; //multipole ak
               }
            }
          }//end if(l==n)
          r=l;
          k = node1->level;

          l=preorder(node1->first_child,l+1,real(node1->loc)-a/(2^(k+1)),imag(node1->loc)-a/(2^(k+1)),node1->loc);
          l=preorder(node1->second_child,l+1,real(node1->loc) + a/(2^(k+1)),imag(node1->loc) + a/(2^(k+1)),node1->loc);
          l=preorder(node1->third_child,l+1,real(node1->loc) - a/(2^(k+1)),imag(node1->loc) - a/(2^(k+1)),node1->loc);
          l=preorder(node1->fourth_child,l+1,real(node1->loc) + a/(2^(k+1)),imag(node1->loc) + a/(2^(k+1)),node1->loc);
          //translation of multipole expansion
          //step 2
          if(r==n){//copy ak to bk
            for(d=1;d<=p;d++){
               node1->bk[d] = node1->ak[d];
            }
          }
          else if(r<n){
               for(f=1;f<=p;f++){
                  for(k=1;k<=r;k++){
                      sec_term1 = sec_term1 + node1->first_child->bk[k]*pow(node1->first_child->Zo,r-k)*binom_coeff(r-1,k-1);
                      sec_term2 = sec_term2 + node1->second_child->bk[k]*pow(node1->second_child->Zo,r-k)*binom_coeff(r-1,k-1);
                      sec_term3 = sec_term3 + node1->third_child->bk[k]*pow(node1->third_child->Zo,r-k)*binom_coeff(r-1,k-1);
                      sec_term4 = sec_term4 + node1->fourth_child->bk[k]*pow(node1->fourth_child->Zo,r-k)*binom_coeff(r-1,k-1);
                  }
                  node1->bk[f] = (node1->first_child->bk[0]*pow(node1->first_child->Zo,r)/(complex<double>)(r)+sec_term1)+
                                 (node1->second_child->bk[0]*pow(node1->second_child->Zo,r)/(complex<double>)(r)+sec_term2)+
                                 (node1->third_child->bk[0]*pow(node1->third_child->Zo,r)/(complex<double>)(r)+sec_term3)+
                                 (node1->fourth_child->bk[0]*pow(node1->fourth_child->Zo,r)/(complex<double>)(r)+sec_term4);
               }          
          }
          return(r-1);
     }
}




int sum_gp(int b){
    int sum = 0,i;
    for(i=0;i<b;i++)
       sum = sum + pow(4,i);
    return(sum);
}
    
//step 3
//level order not recursive but iterative
void levelorder(node* root){
    int j = 0,sum = 0,l = 0,i = 0,b = 1; 
    //node* array[350] = {NULL};//level 4(341)
    if(root == NULL) return;
    Queue_CA Q;
    Q.Enqueue(root);
    while(!Q.IsEmpty()){
         node* current = Q.Front();
         Q.Dequeue();
         j++;
         if(current->first_child != NULL) Q.Enqueue(current->first_child);
         if(current->second_child != NULL) Q.Enqueue(current->second_child);
         if(current->third_child != NULL) Q.Enqueue(current->third_child);
         if(current->fourth_child != NULL) Q.Enqueue(current->fourth_child);
         
         // l = l(j)
         if(j == sum_gp(0))
            l = 1;

         while(1){
            if(j >= sum_gp(b) && j < sum_gp(b+1))
               l = b + 1;
            else
               b = b + 1;
         }
                 
        ////////////
         for(i=0;i<l;i++){
            sum=sum + (4^i);
         }

         
         if(sum == j){// Q has all the nodes of level l
           //do operation on current node
           if(j>=(1)){//then only do//to make sure we are doing from level 2
             //at this point in Q all nodes at level 2 are avialable
             step3func(Q,l); //has to be evaluated from level 1
           }
           if(j>=1+4){//has to be evaluated from level 2
             //if(l<n){//according to the algorithm
             step3_il(Q,l); //Q conatins the nodes which are @ level l and whose il has to be found        
             //il now available
             //now calculate local expansions due to il
             step3_lexp(Q,l);
             //}
             //lexp nw available
           }
           //else{
             
           //end 
         }//end if

         
    }
}
          
//finding nearest neighbours   
void step3func(Queue_CA Q,int l){
    int i,j;
    for(i=0;i<Q.Size();i++){
       int index = (Q.front+i) % MAX_SIZE; 
       for(j=0;j<Q.Size();j++){
          if(i!=j){
            if(abs((Q.A[index])->loc - (Q.A[j])->loc) < 1.4143*a/(2^l)){
               (Q.A[index]->nn).push(Q.A[j]);
            }
           }
       }
    }  //end- for(i
}


//finding interaction list
void step3_il(Queue_CA Q,int l){
 if(l==0 || l==1) return;
 while(!Q.IsEmpty()){
     queue<node*> q_pnn = (Q.Front())->parent->nn;
   while(!q_pnn.empty()){
     if(abs(((q_pnn.front())->first_child->loc)-((Q.Front())->loc)) > 1.4143*a/(2^l))
       ((Q.Front())->il).push((q_pnn.front())->first_child);
     if(abs(((q_pnn.front())->second_child->loc)-((Q.Front())->loc)) > 1.4143*a/(2^l))
       ((Q.Front())->il).push((q_pnn.front())->first_child);
     if(abs(((q_pnn.front())->third_child->loc)-((Q.Front())->loc)) > 1.4143*a/(2^l))
       ((Q.Front())->il).push((q_pnn.front())->first_child);
     if(abs(((q_pnn.front())->fourth_child->loc)-((Q.Front())->loc)) > 1.4143*a/(2^l))
       ((Q.Front())->il).push((q_pnn.front())->first_child);
     q_pnn.pop();
   }//2nd while;Q.front il finding is done
   Q.Dequeue();
 }//1st while  
 return;   
}


//evaluating local expansions
void step3_lexp(Queue_CA Q,int l){
  complex<double> sec_term1 = 0+0I,sec_term2 = 0+0I,one = -1+0I;
  int k,t;
  while(!Q.IsEmpty()){
     node* current = Q.Front();//for 1 such node
     queue<node*> current_ilQ = current->il;
     complex<double> Zo;
     while(!(current_ilQ).empty()){
          node* c_il_node = current_ilQ.front();
          Zo = (current_ilQ.front())->loc - current->loc;
          for(k=1;k<=p;k++)
             sec_term1 = sec_term1 + c_il_node->ak[k]*pow(one,k)/pow(Zo,k);
          current->local_ex_bk[0] = current->local_ex_bk[0] + c_il_node->ak[0]*log(-(Zo)) + sec_term1;
          for(k=1;k<=p;k++){
             for(t=1;t<=p;t++)
                sec_term2 = sec_term2 + c_il_node->ak[t]*pow(one,t)*binom_coeff(k+t-1,t-1)/pow(Zo,t);
             current->local_ex_bk[k] = current->local_ex_bk[k] - c_il_node->ak[0]/((complex<double>)k*pow(Zo,k)) + sec_term2/pow(Zo,k);
          }
          current_ilQ.pop();          
     }//end of 2nd while
     //add parent's local expansion-translated to local_ex_bk
     for(k=0;k<=p;k++)
        for(t=k;t<=p;t++)//for a single k
           current->local_ex_bk[k] = current->local_ex_bk[k] + (current->parent->local_ex_bk[t])*binom_coeff(t,k)*pow((-(current->loc - current->parent->loc)),t-k);

  Q.Dequeue();
  }//end of 1st while
}


complex<double> binom_coeff(int n,int r){
      int i,num = 1;
      for(i=n;i>=n-r+1;i--)
         num = num*i;
      return((complex<double>)(num/fact(r)));
}

int fact(int r){
   int i,f = 1;
   for(i=r;i>=1;i--)
      f = f*i;
   return(f);
}


//now find the potential at any point by using the local expansion formula at the finest level box in which the the point of interest resides
//so find in which box the point of interest resides
