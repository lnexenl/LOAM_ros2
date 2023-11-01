/*
 * This file is modified basing on the ikd-Tree project (https://github.com/hku-mars/ikd-Tree)
 * The original file is released under the GPL v2 License, see project LICENSE file for details
 */

#ifndef IKDTree_H
#define IKDTree_H
#include <cstdio>
#include <queue>
#include <pthread.h>
#include <chrono>
#include <ctime>
#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <memory>
#include <pcl/point_types.h>

#define EPSS 1e-6
#define MinUnbalancedTreeSize 10
#define MultiThreadRebuildPointNum 1500
#define DoDownsample true
#define ForceRebuildPercentage 0.2
#define QSize 1000000

struct BoxPointType{
    float vertex_min[3];
    float vertex_max[3];
};

enum KDTreeOpts {AddPoint, DeletePoint, DeleteBox, AddBox, DownsampleDelete, PushDown};

enum DeletePointStorageOpts {NotRecord, DeletePointsRecord, MultiThreadRecord};

template <typename T>
class Queue{
    private:
        int head = 0,tail = 0, counter = 0;
        T q[QSize];
        bool is_empty = true;
    public:
        void pop();
        T front();
        T back();
        void clear();
        void push(T op);
        bool empty();
        int size();
};



template<typename PointT>
class KDTree{
public:
    using PointVector = std::vector<PointT, Eigen::aligned_allocator<PointT>>;
    using Ptr = std::shared_ptr<KDTree<PointT>>;
    struct KDTreeNode{
        PointT point;
        uint8_t division_axis;  
        int TreeSize = 1;
        int InvalidPointNum = 0;
        int down_del_num = 0;
        bool point_deleted = false;
        bool tree_deleted = false; 
        bool point_downsample_deleted = false;
        bool tree_downsample_deleted = false;
        bool need_push_down_to_left = false;
        bool need_push_down_to_right = false;
        bool working_flag = false;
        float radius_sq;
        pthread_mutex_t push_down_mutex_lock;
        float node_range_x[2], node_range_y[2], node_range_z[2];   
        KDTreeNode *left_son_ptr = nullptr;
        KDTreeNode *right_son_ptr = nullptr;
        KDTreeNode *father_ptr = nullptr;
        // For paper data record
        float alpha_del;
        float alpha_bal;
    };

    struct OperationLoggerType{
        PointT point;
        BoxPointType boxpoint;
        bool tree_deleted, tree_downsample_deleted;
        KDTreeOpts op;
    };

    struct HeapPoint{
        PointT point;
        float dist = 0.0;
        HeapPoint (PointT p = PointT(), float d = INFINITY){
            this->point = p;
            this->dist = d;
        };
        bool operator < (const HeapPoint &a)const{
            if (fabs(dist - a.dist) < 1e-10) return point.x < a.point.x;
            else return dist < a.dist;
        }    
    };

    class Heap{
        public:
            Heap(int max_capacity = 100){
                cap = max_capacity;
                heap = new HeapPoint[max_capacity];
                heap_size = 0;
            }

            ~Heap(){ delete[] heap;}

            void pop(){
                if (heap_size == 0) return;
                heap[0] = heap[heap_size-1];
                heap_size--;
                MoveDown(0);
                return;
            }

            HeapPoint top(){ return heap[0];}

            void push(HeapPoint point){
                if (heap_size >= cap) return;
                heap[heap_size] = point;
                floatUp(heap_size);
                heap_size++;
                return;
            }

            int size(){ return heap_size;}

            void clear(){ heap_size = 0;}
        private:
            int heap_size = 0;
            int cap = 0;        
            HeapPoint * heap;
            void MoveDown(int heap_index){
                int l = heap_index * 2 + 1;
                HeapPoint tmp = heap[heap_index];
                while (l < heap_size){
                    if (l + 1 < heap_size && heap[l] < heap[l+1]) l++;
                    if (tmp < heap[l]){
                        heap[heap_index] = heap[l];
                        heap_index = l;
                        l = heap_index * 2 + 1;
                    } else break;
                }
                heap[heap_index] = tmp;
                return;
            }
            
            void floatUp(int heap_index){
                int ancestor = (heap_index-1)/2;
                HeapPoint tmp = heap[heap_index];
                while (heap_index > 0){
                    if (heap[ancestor] < tmp){
                        heap[heap_index] = heap[ancestor];
                        heap_index = ancestor;
                        ancestor = (heap_index-1)/2;
                    } else break;
                }
                heap[heap_index] = tmp;
                return;
            }

    };    

private:
    // Multi-thread Tree rebuild
    bool terminationFlag = false;
    bool rebuildFlag = false;
    pthread_t rebuildThread{};
    pthread_mutex_t termination_flag_mutex_lock{}, rebuild_ptr_mutex_lock{}, working_flag_mutex{}, search_flag_mutex{};
    pthread_mutex_t rebuild_logger_mutex_lock{}, points_deleted_rebuild_mutex_lock{};
    // queue<OperationLoggerType> rebuildLogger;
    Queue<OperationLoggerType> rebuildLogger;    
    PointVector rebuildPCLStorage;
    KDTreeNode ** Rebuild_Ptr = nullptr;
    int search_mutex_counter = 0;
    static void * multi_thread_ptr(void *arg);
    void multi_thread_rebuild();
    void start_thread();
    void stop_thread();
    void run_operation(KDTreeNode ** root, OperationLoggerType operation);
    // KD Tree Functions and augmented variables
    int tmpTreeSize = 0, tmpValidNum = 0;
    float alpha_bal_tmp = 0.5, alpha_del_tmp = 0.0;
    float deleteParam = 0.5f;
    float balanceParam = 0.7f;
    float downsampleSize = 0.2f;
    bool delStorageDisabled = false;
    KDTreeNode * STATIC_ROOT_NODE = nullptr;
    PointVector deletedPoints;
    PointVector downsampleStorage;
    PointVector deletedPointsMultithread;
    void init_tree_node(KDTreeNode * root);
    void test_lock_states(KDTreeNode *root);
    void buildTree(KDTreeNode ** root, int l, int r, PointVector & Storage);
    void rebuild(KDTreeNode ** root);
    int range_delete(KDTreeNode ** root, BoxPointType boxpoint, bool allow_rebuild, bool is_downsample);
    void point_delete(KDTreeNode ** root, PointT point, bool allow_rebuild);
    void point_add(KDTreeNode ** root, PointT point, bool allow_rebuild, int father_axis);
    void range_add(KDTreeNode ** root, BoxPointType boxpoint, bool allow_rebuild);
    void search(KDTreeNode * root, int k_nearest, PointT point, Heap &q, double max_dist);//priority_queue<HeapPoint>
    void range_search(KDTreeNode *root, BoxPointType boxpoint, PointVector &Storage);
    void radius_search(KDTreeNode *root, PointT point, float radius, PointVector &Storage);
    bool criterion_check(KDTreeNode * root);
    void push_down(KDTreeNode * root);
    void Update(KDTreeNode * root); 
    void delete_tree_nodes(KDTreeNode ** root);
    void downsample(KDTreeNode ** root);
    bool is_same_point(PointT a, PointT b);
    float calc_dist(PointT a, PointT b);
    float calc_box_dist(KDTreeNode * node, PointT point);    
    static bool point_cmp_x(PointT a, PointT b); 
    static bool point_cmp_y(PointT a, PointT b); 
    static bool point_cmp_z(PointT a, PointT b); 

public:
    KDTree(float delete_param = 0.5, float balance_param = 0.6 , float box_length = 0.2);
    ~KDTree();
    void set_delete_param(float delete_param);
    void Set_balance_param(float balance_param);
    void set_downsample_param(float box_length);
    void initialize_KDTree(float _deleteParam = 0.5, float _balanceParam = 0.7, float box_length = 0.2);
    int size();
    int validnum();
    void root_alpha(float &alpha_bal, float &alpha_del);
    void build(PointVector point_cloud);
    void nearest_search(PointT point, int k_nearest, PointVector &Nearest_Points, std::vector<float> & Point_Distance, double max_dist = INFINITY);
    void box_search(const BoxPointType &Box_of_Point, PointVector &Storage);
    void radius_search(PointT point, const float radius, PointVector &Storage);
    int add_points(PointVector & PointToAdd, bool downsample_on);
    void add_point_boxes(std::vector<BoxPointType> & BoxPoints);
    void delete_points(PointVector & PointToDel);
    int delete_point_boxes(std::vector<BoxPointType> & BoxPoints);
    void flatten(KDTreeNode * root, PointVector &Storage, DeletePointStorageOpts storage_type);
    void acquire_removed_points(PointVector & removed_points);
    BoxPointType tree_range();
    PointVector PCLStorage;
    KDTreeNode * rootNode = nullptr;
    int maxQueueSize = 0;
};
#endif
// template <typename PointT>
// PointT KDTree<PointT>::zeroP = PointT(0,0,0);
