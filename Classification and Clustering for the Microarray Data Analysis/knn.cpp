#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;

//K-nearest neigbor classification을 수행하는 함수이다.
//positive train set, negative train set, test set, K, p를 인자로 갖고 test set의 각 row에 대한 prediction 결과를 담은 vector를 반환한다.
vector<string> knn(vector<vector<double>> p_train,  vector<vector<double>> n_train, vector<vector<double>>test , double K, double p){
    vector<string> predictions;
    //test set의 모든 gene마다 모든 train set의 gene과 유전자 발현 데이터의 euclidean distance를 계산한다.
    for(int i=0; i<test.size();i++){
        vector<double> distances;
        double positive = 0;
        double negative = 0;
        for(int j=0; j<p_train.size();j++){
            double square_distance = 0;
            for(int k=0; k<p_train[j].size();k++){
                square_distance += pow(test[i][k]-p_train[j][k], 2);
            }
            distances.push_back(sqrt(square_distance));
        }
        for(int j=0; j<n_train.size();j++){
            double square_distance = 0;
            for(int k=0; k<n_train[j].size();k++){
                square_distance += pow(test[i][k]-n_train[j][k],2);
            }
            distances.push_back(sqrt(square_distance));
        }
        //가장 euclidean distance가 작은 K개의 gene을 보고 ribosomal(positive)인 것과 non-ribosomal(negative)인 것의 개수를 센다.
        for(int j=0; j<K;j++){
            int min_idx = min_element(distances.begin(), distances.end())-distances.begin();
            if(min_idx < p_train.size()) positive++;
            else negative++;
            distances.erase(distances.begin()+ min_idx);
        }
        //p% 이상이 ribosomal gene이면 해당 test data를 ribosomal gene이라고 분류한다.
        //p% 미만이 ribosomal gene이면 해당 test data를 non-ribosomal gene이라고 분류한다.
        if(positive/K >= p) predictions.push_back("ribosomal");
        else predictions.push_back("non-ribosomal");
    }
    //test set의 각 gene에 대한 prediction 결과를 담은 predictions vector를 반환한다.
    return predictions;
}

//double type의 data가 저장된 1차원 vector를 인자로 받아 해당 vector의 모든 원소의 평균을 계산해서 그 평균을 반환하는 함수이다.
double mean(vector<double>result){
    double n = 0;
    for(int i = 0; i < result.size(); i++){
        n += result[i];
    }
    return n/result.size();
}
int main(){

    //ribosomal gene들의 발현 데이터를 담은 행 121개, 열 79개인 vector를 만든다.
    vector<vector<double>>ribo(121,vector<double>(79));
    //non-ribosomal gene들의 발현 데이터를 담은 행 2346개, 열 79개인 vector를 만든다.
    vector<vector<double>>non_ribo(2346,vector<double>(79));
    
    //ribo-data.txt file에서 데이터를 읽어 ribo vector에 저장한다.
    ifstream f_ribo;
    f_ribo.open("ribo-data.txt");
    vector<vector <string>> str_ribo;
    string str;
    while(getline(f_ribo, str)){
        stringstream sstream(str);
        vector<string> gene;
        string sample;
        while(getline(sstream,sample,'\t')){
            gene.push_back(sample);
        }
        str_ribo.push_back(gene);
    }
    //file에서 읽은 데이터는 string type이므로 double type으로 바꿔 ribo vector에 저장한다.
    for (int i=0; i<121; i++){
        for(int j=0; j<79;j++){
            ribo[i][j] = stod(str_ribo[i][j]);
        }
    }
    
    //nonribo-data.txt file에서 데이터를 읽어 non-ribo vector에 저장한다.
    ifstream f_non_ribo;
    f_non_ribo.open("nonribo-data.txt");
    vector<vector <string>> str_non_ribo;
    while(getline(f_non_ribo, str)){
        stringstream sstream(str);
        vector<string> gene;
        string sample;
        while(getline(sstream,sample,'\t')){
            gene.push_back(sample);
        }
        str_non_ribo.push_back(gene);
    }
    //file에서 읽은 데이터는 string type이므로 double type으로 바꿔 non_ribo vector에 저장한다.
    for (int i=0; i<2346; i++){
        for(int j=0; j<79;j++){
            non_ribo[i][j] = stod(str_non_ribo[i][j]);
        }
    }

    //사용자로부터 K와 p를 입력받는다.
    double K;
    double p;
    cout << "K 와 p를 입력하세요." << endl;
    cout << "K : ";
    cin >> K;
    cout << "p : ";
    cin >> p;
    
    double tp = 0; //true positive
    double tn = 0; //true negative
    double fp = 0; //false positive
    double fn = 0; //false negative

    //non-ribosomal gene들 중에 ribosomal gene으로 잘못 분류되는 유전자 번호가 담긴 vector이다.
    vector<int> fp_gene_numbers;
    
    //6-fold cross-validation을 수행한다.
    for(int i=0; i<6;i++){
        vector<vector<double>>p_train = ribo; //positive train set
        vector<vector<double>>n_train = non_ribo; //negative train set
        vector<vector<double>>test; //test set

        //(i+1)번째 fold가 test set으로 쓰이고, 나머지 5개의 fold들이 train set으로 쓰인다.
        //positive set에서 6번째 fold만 gene의 개수가 121개이고 나머지 fold들은 120개이므로 if문으로 6번째 fold인 경우는 따로 처리해준다.
        //negative set은 하나의 fold에 있는 gene의 개수가 391개이다.
        if(i!=5){
            //(i+1)번째 fold를 제거함으로써 positive train set을 만든다.
            p_train.erase(p_train.begin()+i*20, p_train.begin()+(i+1)*20);
            //(i+1)번째 fold의 data를 test vector에 추가한다.
            for(int j = 0; j < 20;j++){
                test.push_back(ribo[i*20+j]);
            }
        }
        else{
            p_train.erase(p_train.begin()+5*20, p_train.begin()+(5+1)*20+1);
            for(int j = 0; j < 21;j++){
                test.push_back(ribo[i*20+j]);
            }
        }
        //(i+1)번째 fold를 제거함으로써 negative train set을 만든다.
        n_train.erase(n_train.begin()+i*391, n_train.begin()+(i+1)*391);
        //(i+1)번째 fold의 data를 test vector에 추가한다.
        for(int j = 0; j < 391;j++){
                test.push_back(non_ribo[i*391+j]);
        }
        //만들어진 positive train set, negative train set, test set, 그리고 사용자가 지정한 K와 p value를 인자로 넣고 
        //K-nearest neigbor classification을 수행하고 예측한 classification 결과를 predictions vector에 저장한다.
        vector<string> predictions = knn(p_train, n_train, test, K, p);

        //모델이 prediction한 결과에 따라 true positive, false negative, true negative, false positive의 개수를 구한다.
        for(int j=0; j<predictions.size();j++){
            int ribo_size = test.size()-391;
            if(j < ribo_size){
                if(predictions[j]=="ribosomal") tp++;
                else fn++;
            }
            else{
                if(predictions[j]=="non-ribosomal") tn++;
                else {
                    fp++;
                    //non-ribosomal gene인데 ribosomal gene으로 잘못 분류된 유전자 번호를 fp_gene_numbers vector에 저장한다.
                    fp_gene_numbers.push_back(391*i + (j-ribo_size)+1);
                }
            }
        }
    }

    //6개의 test 결과로 나온 tp, fn, tn, fp를 합하여 최종 sensitivity, specificity, accuracy를 구한다.
    double sensitivity = tp/(tp+fn);
    double specificity = tn/(tn+fp);
    double accuracy = (tp+tn)/(tp+tn+fp+fn);

    //만든 output 문자열들을 terminal screen에서 보여준다.
    cout << endl;
    cout << "k: " << fixed << setprecision(0) << K << endl;
    cout << "p: " << fixed << setprecision(2) << p << endl;
    cout << "sensitivity: " <<  fixed << setprecision(2) << sensitivity << endl;
    cout << "specificity: "<< fixed << setprecision(2) <<specificity << endl;
    cout << "accuracy: " << fixed << setprecision(2) <<accuracy << endl;

    //만든 output 문자열들을 "knn.out" file에 저장한다.
    ofstream f_result;
    f_result.open("knn.out");
    if(f_result.is_open()){
        f_result << "k: " << fixed << setprecision(0) << K << endl;
        f_result << "p: " << fixed << setprecision(2) << p << endl;
        f_result << "sensitivity: " << fixed << setprecision(2) << sensitivity << endl;
        f_result << "specificity: "<< fixed << setprecision(2) <<specificity << endl;
        f_result << "accuracy: " << fixed << setprecision(2) <<accuracy << endl;
    }

}