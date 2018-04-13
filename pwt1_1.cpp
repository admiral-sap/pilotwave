#ifndef _CLESS_
#define _CLESS_
#include <iostream>

class Board{
	double* board;
	public:
	
	Board();
	~Board();
	Board(const Board&);
	
	double& operator() (int,int);
	const double& operator() (int,int) const;
	double& operator[](int i);
	const double& operator[](int i) const;
  	
	Board& operator=(const Board&);
	
	void norm();
	void nilNorm();
 
};

const Board operator+(const Board&, const Board&);
const Board operator-(const Board&, const Board&);
const Board operator*(const double& x, const Board& A);
const Board operator*(const Board& A, const double& x);
std::ostream& operator<<(std::ostream& output, const Board& A);
void launchGame();

class Game {
	Board* board; 		//array of all players boards
	int player;			//currently active player
	int playernumber;	//number of players
	bool* alive;			//array of players' status
public:

	Game();									//create two player game
	~Game();								//delete game
	Game(int number, int*pos);		//create gane with 'number' players and xth player's starting position (pos[2*x],pos[2*x+1])

	void readMove();						//read in move from player
	void print();
	void makeMove(int i, int j, int sgn_1, int sgn_2, int k, int mode, int observe);
	void executeMove(Board move);
	void norm();
	bool isFinished();
};

#endif
#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <random>
#include <cmath>

const int BOARDSZ_1=6;
const int BOARDSZ_2=6;
const std::string PLAYERBRACKET_L[8] ={"(","[","{","|","!","'",")","*"};
const std::string PLAYERBRACKET_R[8] ={")","]","}","|","!","'","(","*"};
const double CAREFUL_ACCURACY=0.7;
double violent(int k, int l){
	if(k==l)
		return 5.;
	else
		return (1./(fabs(k-l)));
}

Board::Board(){
	board=new double[BOARDSZ_1*BOARDSZ_2];
	for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++){
		board[i]=0.;
	}
}

Board::~Board(){
	delete[] board;
}

Board::Board(const Board& rhs){
	board=new double[BOARDSZ_1*BOARDSZ_2];
	for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++){
		board[i]=rhs.board[i];
	}
}

double& Board::operator()(int i, int j){
   assert(i<BOARDSZ_1&&j<BOARDSZ_2);
	return board[i+BOARDSZ_1*j];
} 

const double& Board::operator()(int i, int j) const{
 assert(i<BOARDSZ_1&&j<BOARDSZ_2);
	return board[i+BOARDSZ_1*j];
}

double& Board::operator[](int i){
  assert(i<BOARDSZ_1*BOARDSZ_2);
	return board[i];
}

const double& Board::operator[](int i) const{
  assert(i<BOARDSZ_1*BOARDSZ_2);
	return board[i];
}

Board& Board::operator=(const Board& rhs) {
	if ( this != &rhs) {
		for (int ell=0; ell<BOARDSZ_1*BOARDSZ_2; ++ell) {
			board[ell] = rhs[ell];
		}
	}
	return *this;
}

const Board operator+(const Board& lhs, const Board& rhs){
	Board out;
	for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++){
		out[i]=lhs[i]+rhs[i];
	}
	return out;
}

const Board operator*(const Board& A, const double& x){
	Board out;
	for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++){
		out[i]=A[i]*x;
	}
	return out;
}

const Board operator*(const double& x, const Board& A){
	Board out;
	for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++){
		out[i]=A[i]*x;
	}
	return out;
}

const Board operator-(const Board& lhs, const Board& rhs){
	Board out;
	for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++){
		out[i]=lhs[i]-rhs[i];
	}
	return out;
}

void Board::norm(){
  double sum=0.;
  for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++)
    sum+=board[i];
  for(int i=0;i<BOARDSZ_1*BOARDSZ_2;i++)
    board[i]=board[i]/sum;
}

std::ostream& operator<<(std::ostream& output, const Board& A) {
	std::cout<< std::fixed << std::setprecision(2);
	for(int l=0; l<BOARDSZ_2;l++){
			output<<"_______";
		}
	output << "|\n";
	for (int j=0; j<BOARDSZ_1; j++) {
		for (int k=0; k<BOARDSZ_2; k++) {
			output << " | " << A(j,k);
		}
		output << "|\n";
		for(int l=0; l<BOARDSZ_2;l++){
			output<<"_______";
		}
		output<<"|\n";
	}
	output<<"\n\n";
	return output;
}



Game::Game(){
  player=0;
  playernumber=2;
  board=new Board[2];
  alive=new bool[2];
  alive[0]=true;
  alive[1]=true;
  board[0](0,0)=1.;
  board[1](BOARDSZ_1-1,BOARDSZ_2-1)=1.; 
}

Game::~Game(){
  delete[] board;
}

Game::Game(int number, int* pos){
	player=0;
	playernumber=number;
	board=new Board[number];
	for(int i=0;i<playernumber;i++){
		board[i](pos[2*i],pos[2*i+1])=1.;
	}
	alive=new bool[number];
	for(int i=0;i<playernumber;i++){
    alive[i]=true;
 }
}

void Game::readMove(){

  if(alive[player]==true){
	 print();
	int i;
	int j;
	std::string direction;
	int k;
	int sgn_1;
	int sgn_2;
	int mode;
	int observe;
	bool valid=false;
	while(valid==false){ //read in valid move
		valid=true;
		std::cout<<std::endl<<std::endl<<"Player "<<PLAYERBRACKET_L[player]<<player<<PLAYERBRACKET_R[player]<<". Please enter your move."<<std::endl<<"Coordinates of starting square(row,column)"<<std::endl;
		std::cin>>i>>j;
		std::cout<<std::endl<<"Direction (N/NE/E/SE/S/SW/W/NW): ";
		std::cin>>direction;
		std::cout<<std::endl<<"Distance: ";
		std::cin>>k;
		std::cout<<std::endl<<"Mode of execution (1=carefully, 2=violently): ";
		std::cin>>mode;
		std::cout<<std::endl<<"Mode of observation (0=none, 1=total): ";
		std::cin>>observe;
		if (direction=="N"){ 	//translate direction
			sgn_1=-1;
			sgn_2=0;
		}else if(direction=="NE"){
			sgn_1=-1;
			sgn_2=1;
		}else if(direction=="E"){
			sgn_1=0;
			sgn_2=1;
		}else if(direction=="SE"){
			sgn_1=1;
			sgn_2=1;
		}else if(direction=="S"){
			sgn_1=1;
			sgn_2=0;
		}else if(direction=="SW"){
			sgn_1=1;
			sgn_2=-1;
		}else if(direction=="W"){
			sgn_1=0;
			sgn_2=-1;
		}else if(direction=="NW"){
			sgn_1=-1;
			sgn_2=-1;
		}else{
			valid=false;
		}
		if (k<=0){					//assert validity of move
			valid=false;
			std::cout<<std::endl<<"invalid distance";
		} if (i+(sgn_1*k)>=BOARDSZ_1||i+(sgn_1*k)<0){
			valid=false;
			std::cout<<std::endl<<"move exceeds board";
		} if (j+(sgn_2*k)>=BOARDSZ_2||j+(sgn_2*k)<0){
			valid=false;
			std::cout<<std::endl<<"move exceeds board";
		}if ( i<0 || i>=BOARDSZ_1 || j<0 || j>=BOARDSZ_2){
			valid=false;
			std::cout<<std::endl<<"impossible starting sqaure";
		}else if (board[player](i,j)==0.){
			valid=false;
			std::cout<<std::endl<<"invalid starting sqaure";
		} if(mode!=1&&mode!=2){
			valid=false;
			std::cout<<std::endl<<"invalid mode of execution";
		} if(observe!=0&&observe!=1){
			valid=false;
			std::cout<<std::endl<<"invalid mode of observation";
		}
	}
	makeMove(i,j,sgn_1,sgn_2,k,mode,observe);
	}
	if(isFinished()){
		std::cout<<"Player "<<PLAYERBRACKET_L[player]<<player<<PLAYERBRACKET_R[player]<<" won the game\n";
		getchar();
		std::cout<<"Press any key to end ...";
		getchar();
	}else {
		player=(player+1)%playernumber;
		readMove();
	}
}

void Game::makeMove(int i, int j, int sgn_1, int sgn_2, int k, int mode, int observe){
	if(board[player](i,j)>0){
	
	unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	Board dist_pos;
	Board dist_neg;
	Board move;
	double val;
	bool roll_success=false;
	int enemy=-1;
	
	for(int p=0;p<playernumber;p++){
			if(board[p](i+(k*sgn_1),j+(k*sgn_2))>0){
				if(p!=player){
				enemy=p;
				break;
				}
			}
	}
	if (mode==1){
		val=board[player](i,j);
		dist_neg(i,j)=1.;
		dist_pos(i+(k*sgn_1),j+(k*sgn_2))=CAREFUL_ACCURACY;
		dist_pos(i,j)=1.-CAREFUL_ACCURACY;
		if(observe==0){
			move=(dist_pos-dist_neg)*val;
			executeMove(move);
		}
		else if(observe==1){
			if(distribution(generator)<dist_pos(i+(k*sgn_1),j+(k*sgn_2))){
				Board tmp;
				tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
				dist_pos=tmp;
				if(distribution(generator)<val+board[player](i+(k*sgn_1),j+(k*sgn_2))){
					Board tmp;
					tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
					board[player]=tmp;
					if(enemy!=-1){
						if(distribution(generator)<board[enemy](i+(k*sgn_1),j+(k*sgn_2))){
							Board tmp;
							board[enemy]=tmp;
							alive[enemy]=false;
						}
						else{
							board[enemy](i+(k*sgn_1),j+(k*sgn_2))=0;
						}
					}
					norm();
				}
				else{
					board[player](i+(k*sgn_1),j+(k*sgn_2))=0.;
					move=dist_neg * val * -1.;
					executeMove(move);
					if(enemy!=-1){
						if(distribution(generator)<board[enemy](i+(k*sgn_1),j+(k*sgn_2))){
							Board tmp;
							tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
							board[enemy]=tmp;
						}
						else{
							board[enemy](i+(k*sgn_1),j+(k*sgn_2))=0.;
						}
					}
					norm();
				}
			}
			else{
				dist_pos(i+(k*sgn_1),j+(k*sgn_2))=0.;
				dist_pos.norm();
				move=(dist_pos-dist_neg)*val;
				executeMove(move);
			}
			
		}
	}
	else if(mode==2){
		if(sgn_1==0){
			if(i+1<BOARDSZ_1)
				makeMove(i+1, j,sgn_1,sgn_2,k, 0, 0);
			if(i-1>=0)
				makeMove(i-1, j,sgn_1,sgn_2,k, 0, 0);
		}else if(sgn_2==0){
			if(j+1<BOARDSZ_2)
				makeMove(i, j+1,sgn_1,sgn_2,k, 0, 0);
			if(j-1>=0)
				makeMove(i, j-1,sgn_1,sgn_2,k, 0, 0);
		}
		val=board[player](i,j);
		dist_neg(i,j)=1.;
		for(int l=0;l<=2*k;l++){
			if(i+(l*sgn_1)<BOARDSZ_1 && j+(l*sgn_2)<BOARDSZ_2 && i+(l*sgn_1)>=0 && j+(l*sgn_2)>=0)
				dist_pos(i+(l*sgn_1),j+(l*sgn_2))=violent(k,l);
		}
		dist_pos.norm();
		if(observe==0){
			move=(dist_pos-dist_neg)*val;
			executeMove(move);
		}
		else if(observe==1){
			if(distribution(generator)<dist_pos(i+(k*sgn_1),j+(k*sgn_2))){
				Board tmp;
				tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
				dist_pos=tmp;
				if(distribution(generator)<val+board[player](i+(k*sgn_1),j+(k*sgn_2))){
					Board tmp;
					tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
					board[player]=tmp;
					if(enemy!=-1){
						if(distribution(generator)<board[enemy](i+(k*sgn_1),j+(k*sgn_2))){
							Board tmp;
							board[enemy]=tmp;
							alive[enemy]=false;
						}
						else{
							board[enemy](i+(k*sgn_1),j+(k*sgn_2))=0;
						}
					}
					norm();
				}
				else{
					board[player](i+(k*sgn_1),j+(k*sgn_2))=0.;
					move=dist_neg * val * -1.;
					executeMove(move);
					if(enemy!=-1){
						if(distribution(generator)<board[enemy](i+(k*sgn_1),j+(k*sgn_2))){
							Board tmp;
							tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
							board[enemy]=tmp;
						}
						else{
							board[enemy](i+(k*sgn_1),j+(k*sgn_2))=0.;
						}
					}
					norm();
				}
			}
			else{
				dist_pos(i+(k*sgn_1),j+(k*sgn_2))=0.;
				dist_pos.norm();
				move=(dist_pos-dist_neg)*val;
				executeMove(move);
			}
			
		}
	}
	
	else if(mode==0){
		val=board[player](i,j);
		dist_neg(i,j)=1.;
		for(int l=0;l<=2*k;l++){
			if(i+(l*sgn_1)<BOARDSZ_1 && j+(l*sgn_2)<BOARDSZ_2)
				dist_pos(i+(l*sgn_1),j+(l*sgn_2))=violent(k,l);
		}
		dist_pos.norm();
		if(observe==0){
			move=(dist_pos-dist_neg)*val;
			executeMove(move);
		}
		else if(observe==1){
			if(distribution(generator)<dist_pos(i+(k*sgn_1),j+(k*sgn_2))){
				Board tmp;
				tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
				dist_pos=tmp;
				if(distribution(generator)<val+board[player](i+(k*sgn_1),j+(k*sgn_2))){
					Board tmp;
					tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
					board[player]=tmp;
					if(enemy!=-1){
						if(distribution(generator)<board[enemy](i+(k*sgn_1),j+(k*sgn_2))){
							Board tmp;
							board[enemy]=tmp;
							alive[enemy]=false;
						}
						else{
							board[enemy](i+(k*sgn_1),j+(k*sgn_2))=0;
						}
					}
					norm();
				}
				else{
					board[player](i+(k*sgn_1),j+(k*sgn_2))=0.;
					move=dist_neg * val * -1.;
					executeMove(move);
					if(enemy!=-1){
						if(distribution(generator)<board[enemy](i+(k*sgn_1),j+(k*sgn_2))){
							Board tmp;
							tmp(i+(k*sgn_1),j+(k*sgn_2))=1.;
							board[enemy]=tmp;
						}
						else{
							board[enemy](i+(k*sgn_1),j+(k*sgn_2))=0.;
						}
					}
					norm();
				}
			}
			else{
				dist_pos(i+(k*sgn_1),j+(k*sgn_2))=0.;
				dist_pos.norm();
				move=(dist_pos-dist_neg)*val;
				executeMove(move);
			}
			
		}
	}
}
	
	
	
	/*if(isFinished()){
		std::cout<<"Player "<<PLAYERBRACKET_L[player]<<player<<PLAYERBRACKET_R[player]<<" won the game";
	}else if(mode!=0){
		player=(player+1)%playernumber;
		readMove();
	}*/
}

void Game::executeMove(Board move){
	for (int i=0; i<BOARDSZ_1*BOARDSZ_2;i++){
		for(int p=0;p<playernumber;p++){
			if(p==player) {
				
			}
			else if(board[p][i]>0){
				if(board[p][i]>move[i]){
					board[p][i]=board[p][i]-move[i];
					move[i]=0.;
				}else{
					move[i]=move[i]-board[p][i];
					board[p][i]=0.;
				}
			}
		}
	}
	board[player]=board[player]+move;
	norm();
}
void Game::print(){
	std::cout<< std::fixed << std::setprecision(2);
	for(int l=0; l<BOARDSZ_2;l++){
			std::cout<<"_______";
		}
	std::cout << "\n";
	for (int j=0; j<BOARDSZ_1; j++) {
		for (int k=0; k<BOARDSZ_2; k++) {
			std::cout << "|";
			for(int i=0; i<=playernumber;i++){
				if(i==playernumber){
					std::cout<<" "<<"    "<<" ";
				}
				else if(board[i](j,k)>0){
					std::cout<<PLAYERBRACKET_L[i]<<board[i](j,k)<<PLAYERBRACKET_R[i];
					break;
				}
			}
		}
		std::cout << "|\n";
		for(int l=0; l<BOARDSZ_2;l++){
			std::cout<<"_______";
		}
		std::cout<<"|\n";
	}
	std::cout<<"\n\n";
}

void Game::norm(){
	for(int i=0;i<playernumber;i++)
		if (alive[i]==true)
			board[i].norm();
}

bool Game::isFinished(){
	int living=0;
	int winner=-1;
	for(int i=0;i<playernumber;i++){
		if(alive[i]==true){
			living++;
			winner=i;
		}
	}
	if(living<=1){
		return true;
	}
	else{
		return false;
	}
}
void launchGame(){
	bool valid;
	int playernumber;
	std::cout<<std::endl<<std::endl<<"Pilot Wave: The game | Version: 1.1"<<std::endl<<"an Aniline production"<<std::endl<<std::endl<<"Number of players:"; // Enter Version Number here (should be changed to a String Version Number on top)
	std::cin>>playernumber;
	int* pos=new int[2*playernumber];
		for(int i=0;i<playernumber;i++){
			valid=false;
			while (valid==false){
				valid=true;
			std::cout<<"Player "<<PLAYERBRACKET_L[i]<<i<<PLAYERBRACKET_R[i]<<", please enter your starting position (row/column):"<<std::endl;
			std::cin>>pos[2*i]>>pos[(2*i)+1];
			if (pos[2*i]<0 || pos[2*i]>=BOARDSZ_1)
				valid=false;
			if (pos[(2*i)+1]<0 || pos[(2*i)+1]>=BOARDSZ_2)
				valid=false;
			}
		}
		Game newgame(playernumber,pos);
		newgame.readMove();
}

int main(){
	launchGame();
}