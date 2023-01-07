#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <ctime>
using namespace std;

vector<int> knn(vector<vector<double>> dataset, int K, vector<int> center_points){

    int r_size = dataset.size(); //dataset의 행의 개수
    int c_size = dataset[0].size(); //dataset의 열의 개수

    vector<int>result(r_size);
    vector<vector<double>> distances(K, vector<double>(r_size));
    vector<vector<double>> centers(K, vector<double>(c_size));

    //center_points 인자에 따라 K개의 center를 선택한다.
    for(int i = 0; i < K; i++){
        int idx = center_points[i]-1;
        for(int j = 0; j < c_size; j++){
            centers[i][j] = dataset[idx][j];
        }
    }

    //K개의 cluster center와 dataset과의 거리를 계산해서 그 중 가장 가까운 cluster center의 cluster로 dataset의 각각의 data point를 assign한다.
    //각각의 data point가 K개의 cluster로 분할이 된다.
    for(int i = 0; i < r_size; i++){
        double min;
        for(int j = 0; j < K; j++){
            double square_distance = 0;
            for(int k = 0; k < c_size; k++){
                square_distance += pow(dataset[i][k]-centers[j][k],2);
            }
            distances[j][i] = sqrt(square_distance);
            if(j==0) min = distances[j][i];
            if(distances[j][i]<=min){
                result[i] = j+1;
                min = distances[j][i];
            }
        }
    }

    //K개의 cluster center를 다시 계산하고 data point를 업데이트된 cluster center를 기준으로 거리를 다시 계산해서 그 중 가장 가까운 cluster로 data point들을 assign한다.
    //clustering 결과가 달라지지 않을 때까지 반복한다.
    bool is_change = true; 
    while(is_change){
        is_change = false;
        //K개의 cluster center를 다시 계산한다.
        vector<vector<double>> next_centers(K, vector<double>(c_size, 0));
        vector<double> num(K, 0);
   
        for(int i = 0; i < r_size;i++){
            int idx = result[i]-1;
            num[idx]+=1;
            for(int j = 0; j < c_size;j++){
                next_centers[idx][j] += dataset[i][j];
            }
        }

        for(int i = 0; i < K;i++){
            for(int j = 0; j < c_size;j++){
                next_centers[i][j]=next_centers[i][j]/num[i];
                centers[i][j] = next_centers[i][j];
            }
        }
        //업데이트된 cluster center를 기준으로 거리를 다시 계산해서 그 중 가장 가까운 cluster로 data point들을 assign한다.
        //하나라도 clustering 결과가 달라진 data point가 있으면 while문을 반복한다.
        vector<int>update_result(r_size);
        for(int i = 0; i < r_size; i++){
            double min = 1000000;
            for(int j = 0; j < K; j++){
                double square_distance = 0;
                for(int k = 0; k < c_size; k++){
                    square_distance += pow(dataset[i][k]-centers[j][k],2);
                }
                distances[j][i] = sqrt(square_distance);
                if(distances[j][i]<min){
                    update_result[i] = j+1;
                    min = distances[j][i];
                }
            }
            if(update_result[i] != result[i]){
                result[i] = update_result[i];
                is_change = true;
            }
        }
    }

    //while문을 나갔다는 것은 data point들의 cluster가 변경이 없이 수렴했다는 의미이므로 마지막에 저장된 result vector가 K-means clustering 결과이다.
    //해당 result vector를 반환한다.
    return result;
}

int main(){

    //ribosomal gene과 non-ribosomal gene들의 발현 데이터를 담은 행 2467개, 열 79개인 vector를 만든다.
    vector<vector<double>>dataset(2467,vector<double>(79));
    
    //ribo-data.txt file에서 데이터를 읽어 dataset vector에 저장한다.
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
    //file에서 읽은 데이터는 string type이므로 double type으로 바꿔 dataset vector에 저장한다.
    for (int i=0; i<121; i++){
        for(int j=0; j<79;j++){
            dataset[i][j] = stod(str_ribo[i][j]);
        }
    }
    
    //nonribo-data.txt file에서 데이터를 읽어 dataset vector에 저장한다.
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
    //file에서 읽은 데이터는 string type이므로 double type으로 바꿔 dataset vector에 저장한다.
    for (int i=121; i<2467; i++){
        for(int j=0; j<79;j++){
            dataset[i][j] = stod(str_non_ribo[i-121][j]);
        }
    }

    //사용자로부터 cluster의 개수 K를 입력받는다.
    int K;
    cout << "cluster의 개수 K : ";
    cin >> K;

    //random number 생성을 위해 time을 seed로 가지는 srand 함수를 호출한다.
    srand((unsigned int)time(NULL));

    //사용자로부터 cluster의 center들을 입력받는다.
    vector<int> center_points(K);
    int c;
    cout << "starting center로 지정할 data point 번호를 입력하시오 (0 : random data point) : ";
    cin >> c;
    
    for(int i = 0; i < K ; i++){
        center_points[i] = c;
    }
    //center_points가 0일 때, random data point를 선택한다.
    if(center_points[0]==0){
        center_points[0] = rand()%121+1; //1부터 121 사이 정수 중에 random하게 뽑는다.
    }
    if(center_points[1]==0){
        center_points[1] = rand()%2346+122; //122부터 2467 사이 정수 중에 random하게 뽑는다.
    }
    else center_points[1] += 121;

    //dataset과 cluster의 개수 K와 center point들을 인자로 넣고 K-means clustering을 수행한 뒤, 그 결과를 result vector에 저장한다.
    vector<int> result = knn(dataset, K, center_points);

    double tp = 0; //true positive
    double tn = 0; //true negative
    double fp = 0; //false positive
    double fn = 0; //false negative

    vector<int> tp_genes; //true positive gene number
    vector<int> fp_genes; //true negative gene number
    vector<int> tn_genes; //false positive gene number
    vector<int> fn_genes; //false negative gene number

    //clustering 결과에 따라 true positive, false negative, true negative, false positive에 해당하는 유전자 번호를 각각의 vector에 저장한다.
    //clustering 결과에 따라 true positive, false negative, true negative, false positive의 개수를 구한다.
    for(int i = 0; i < 2467; i++){
        if(i<121){
            if(result[i]==1){
                tp_genes.push_back(i+1);
                tp+=1;
            }
            else{
                fn_genes.push_back(i+1);
                fn+=1;
            }
        }
        else{
            if(result[i]==2){
                tn_genes.push_back(i+1);
                tn+=1;
            }
            else{
                fp_genes.push_back(i+1);
                fp+=1;
            }
        }
    }
    //ribosomal gene cluster에 있는 ribosomal gene들의 비율 출력
    cout << endl;
    cout << "ribosomal gene cluster에 있는 ribosomal gene들의 비율 : " << tp/(tp+fp) <<endl<<endl; 

    //nonribosomal gene cluster에 있는 ribosomal gene들의 비율 출력
    cout << "non-ribosomal gene cluster에 있는 ribosomal gene들의 비율 : " << fn/(fn+tn) <<endl<<endl; 

    //non-ribosomal gene으로 잘못 clustering된 ribosomal gene number 출력
    cout << "non-ribosomal gene으로 잘못 clustering된 ribosomal gene number : ";
    if(fn_genes.size() == 0) cout<<"None"<<endl;
    else{
        for(int i = 0; i < fn_genes.size(); i++){
            cout << fn_genes[i] << " ";
        }
        cout << endl;
    }
}