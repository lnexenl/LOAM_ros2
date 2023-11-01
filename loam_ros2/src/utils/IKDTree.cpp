/*
 * This file is modified basing on the ikd-Tree project (https://github.com/hku-mars/ikd-Tree)
 * The original file is released under the GPL v2 License, see project LICENSE file for details
 * Description: ikd-Tree: an incremental k-d tree for robotic applications
 * Author: Yixi Cai
 * email: yixicai@connect.hku.hk
 */
#include "utils/IKDTree.h"

using std::min;
using std::max;

template <typename PointT>
KDTree<PointT>::KDTree(float _deleteParam, float _balanceParam, float _boxLength) {
    deleteParam = _deleteParam;
    balanceParam = _balanceParam;
    downsampleSize = _boxLength;
    rebuildLogger.clear();           
    terminationFlag = false;
    start_thread();
}

template <typename PointT>
KDTree<PointT>::~KDTree()
{
    stop_thread();
    delStorageDisabled = true;
    delete_tree_nodes(&rootNode);
    PointVector ().swap(PCLStorage);
    rebuildLogger.clear();           
}

template <typename PointT>
void KDTree<PointT>::set_delete_param(float delete_param){
    deleteParam = delete_param;
}

template <typename PointT>
void KDTree<PointT>::Set_balance_param(float balance_param){
    balanceParam = balance_param;
}

template <typename PointT>
void KDTree<PointT>::set_downsample_param(float downsample_param){
    downsampleSize = downsample_param;
}

template <typename PointT>
void KDTree<PointT>::initialize_KDTree(float deleteParam, float balanceParam, float box_length){
    set_delete_param(deleteParam);
    Set_balance_param(balanceParam);
    set_downsample_param(box_length);
}

template <typename PointT>
void KDTree<PointT>::init_tree_node(KDTreeNode * root){
    root->point.x = 0.0f;
    root->point.y = 0.0f;
    root->point.z = 0.0f;       
    root->node_range_x[0] = 0.0f;
    root->node_range_x[1] = 0.0f;
    root->node_range_y[0] = 0.0f;
    root->node_range_y[1] = 0.0f;    
    root->node_range_z[0] = 0.0f;
    root->node_range_z[1] = 0.0f;     
    root->division_axis = 0;
    root->father_ptr = nullptr;
    root->left_son_ptr = nullptr;
    root->right_son_ptr = nullptr;
    root->TreeSize = 0;
    root->InvalidPointNum = 0;
    root->down_del_num = 0;
    root->point_deleted = false;
    root->tree_deleted = false;
    root->need_push_down_to_left = false;
    root->need_push_down_to_right = false;
    root->point_downsample_deleted = false;
    root->working_flag = false;
    pthread_mutex_init(&(root->push_down_mutex_lock),nullptr);
}   

template <typename PointT>
int KDTree<PointT>::size(){
    int s = 0;
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
        if (rootNode != nullptr) {
            return rootNode->TreeSize;
        } else {
            return 0;
        }
    } else {
        if (!pthread_mutex_trylock(&working_flag_mutex)){
            s = rootNode->TreeSize;
            pthread_mutex_unlock(&working_flag_mutex);
            return s;
        } else {
            return tmpTreeSize;
        }
    }
}

template <typename PointT>
BoxPointType KDTree<PointT>::tree_range(){
    BoxPointType range;
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
        if (rootNode != nullptr) {
            range.vertex_min[0] = rootNode->node_range_x[0];
            range.vertex_min[1] = rootNode->node_range_y[0];
            range.vertex_min[2] = rootNode->node_range_z[0];
            range.vertex_max[0] = rootNode->node_range_x[1];
            range.vertex_max[1] = rootNode->node_range_y[1];
            range.vertex_max[2] = rootNode->node_range_z[1];
        } else {
            memset(&range, 0, sizeof(range));
        }
    } else {
        if (!pthread_mutex_trylock(&working_flag_mutex)){
            range.vertex_min[0] = rootNode->node_range_x[0];
            range.vertex_min[1] = rootNode->node_range_y[0];
            range.vertex_min[2] = rootNode->node_range_z[0];
            range.vertex_max[0] = rootNode->node_range_x[1];
            range.vertex_max[1] = rootNode->node_range_y[1];
            range.vertex_max[2] = rootNode->node_range_z[1];
            pthread_mutex_unlock(&working_flag_mutex);
        } else {
            memset(&range, 0, sizeof(range));
        }
    }
    return range;
}

template <typename PointT>
int KDTree<PointT>::validnum(){
    int s = 0;
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
        if (rootNode != nullptr)
            return (rootNode->TreeSize - rootNode->InvalidPointNum);
        else 
            return 0;
    } else {
        if (!pthread_mutex_trylock(&working_flag_mutex)){
            s = rootNode->TreeSize-rootNode->InvalidPointNum;
            pthread_mutex_unlock(&working_flag_mutex);
            return s;
        } else {
            return -1;
        }
    }
}

template <typename PointT>
void KDTree<PointT>::root_alpha(float &alpha_bal, float &alpha_del){
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
        alpha_bal = rootNode->alpha_bal;
        alpha_del = rootNode->alpha_del;
        return;
    } else {
        if (!pthread_mutex_trylock(&working_flag_mutex)){
            alpha_bal = rootNode->alpha_bal;
            alpha_del = rootNode->alpha_del;
            pthread_mutex_unlock(&working_flag_mutex);
            return;
        } else {
            alpha_bal = alpha_bal_tmp;
            alpha_del = alpha_del_tmp;      
            return;
        }
    }    
}

template <typename PointT>
void KDTree<PointT>::start_thread(){
    pthread_mutex_init(&termination_flag_mutex_lock, NULL);   
    pthread_mutex_init(&rebuild_ptr_mutex_lock, NULL);     
    pthread_mutex_init(&rebuild_logger_mutex_lock, NULL);
    pthread_mutex_init(&points_deleted_rebuild_mutex_lock, NULL); 
    pthread_mutex_init(&working_flag_mutex, NULL);
    pthread_mutex_init(&search_flag_mutex, NULL);
    pthread_create(&rebuildThread, NULL, multi_thread_ptr, (void*) this);
    printf("Multi thread started \n");    
}

template <typename PointT>
void KDTree<PointT>::stop_thread(){
    pthread_mutex_lock(&termination_flag_mutex_lock);
    terminationFlag = true;
    pthread_mutex_unlock(&termination_flag_mutex_lock);
    if (rebuildThread) pthread_join(rebuildThread, NULL);
    pthread_mutex_destroy(&termination_flag_mutex_lock);
    pthread_mutex_destroy(&rebuild_logger_mutex_lock);
    pthread_mutex_destroy(&rebuild_ptr_mutex_lock);
    pthread_mutex_destroy(&points_deleted_rebuild_mutex_lock);
    pthread_mutex_destroy(&working_flag_mutex);
    pthread_mutex_destroy(&search_flag_mutex);     
}

template <typename PointT>
void * KDTree<PointT>::multi_thread_ptr(void * arg){
    KDTree * handle = (KDTree*) arg;
    handle->multi_thread_rebuild();
    return nullptr;    
}

template <typename PointT>
void KDTree<PointT>::multi_thread_rebuild(){
    bool terminated = false;
    KDTreeNode * father_ptr, ** new_node_ptr;
    pthread_mutex_lock(&termination_flag_mutex_lock);
    terminated = terminationFlag;
    pthread_mutex_unlock(&termination_flag_mutex_lock);
    while (!terminated){
        pthread_mutex_lock(&rebuild_ptr_mutex_lock);
        pthread_mutex_lock(&working_flag_mutex);
        if (Rebuild_Ptr != nullptr ){                    
            /* Traverse and copy */
            if (!rebuildLogger.empty()){
                printf("\n\n\n\n\n\n\n\n\n\n\n ERROR!!! \n\n\n\n\n\n\n\n\n");
            }
            rebuildFlag = true;
            if (*Rebuild_Ptr == rootNode) {
                tmpTreeSize = rootNode->TreeSize;
                tmpValidNum = rootNode->TreeSize - rootNode->InvalidPointNum;
                alpha_bal_tmp = rootNode->alpha_bal;
                alpha_del_tmp = rootNode->alpha_del;
            }
            KDTreeNode * old_root_node = (*Rebuild_Ptr);                            
            father_ptr = (*Rebuild_Ptr)->father_ptr;  
            PointVector ().swap(rebuildPCLStorage);
            // Lock search
            pthread_mutex_lock(&search_flag_mutex);
            while (search_mutex_counter != 0){
                pthread_mutex_unlock(&search_flag_mutex);
                usleep(1);             
                pthread_mutex_lock(&search_flag_mutex);
            }
            search_mutex_counter = -1;
            pthread_mutex_unlock(&search_flag_mutex);
            // Lock deleted points cache
            pthread_mutex_lock(&points_deleted_rebuild_mutex_lock);    
            flatten(*Rebuild_Ptr, rebuildPCLStorage, MultiThreadRecord);
            // Unlock deleted points cache
            pthread_mutex_unlock(&points_deleted_rebuild_mutex_lock);
            // Unlock search
            pthread_mutex_lock(&search_flag_mutex);
            search_mutex_counter = 0;
            pthread_mutex_unlock(&search_flag_mutex);              
            pthread_mutex_unlock(&working_flag_mutex);   
            /* rebuild and update missed operations*/
            OperationLoggerType Operation;
            KDTreeNode * new_root_node = nullptr;  
            if (int(rebuildPCLStorage.size()) > 0){
                buildTree(&new_root_node, 0, rebuildPCLStorage.size()-1, rebuildPCLStorage);
                // rebuild has been done. Updates the blocked operations into the new tree
                pthread_mutex_lock(&working_flag_mutex);
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                int tmp_counter = 0;
                while (!rebuildLogger.empty()){
                    Operation = rebuildLogger.front();
                    maxQueueSize = max(maxQueueSize, rebuildLogger.size());
                    rebuildLogger.pop();
                    pthread_mutex_unlock(&rebuild_logger_mutex_lock);                  
                    pthread_mutex_unlock(&working_flag_mutex);
                    run_operation(&new_root_node, Operation);
                    tmp_counter ++;
                    if (tmp_counter % 10 == 0) usleep(1);
                    pthread_mutex_lock(&working_flag_mutex);
                    pthread_mutex_lock(&rebuild_logger_mutex_lock);               
                }   
               pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }  
            /* Replace to original tree*/          
            // pthread_mutex_lock(&working_flag_mutex);
            pthread_mutex_lock(&search_flag_mutex);
            while (search_mutex_counter != 0){
                pthread_mutex_unlock(&search_flag_mutex);
                usleep(1);             
                pthread_mutex_lock(&search_flag_mutex);
            }
            search_mutex_counter = -1;
            pthread_mutex_unlock(&search_flag_mutex);
            if (father_ptr->left_son_ptr == *Rebuild_Ptr) {
                father_ptr->left_son_ptr = new_root_node;
            } else if (father_ptr->right_son_ptr == *Rebuild_Ptr){             
                father_ptr->right_son_ptr = new_root_node;
            } else {
                throw "Error: Father ptr incompatible with current node\n";
            }
            if (new_root_node != nullptr) new_root_node->father_ptr = father_ptr;
            (*Rebuild_Ptr) = new_root_node;
            int valid_old = old_root_node->TreeSize-old_root_node->InvalidPointNum;
            int valid_new = new_root_node->TreeSize-new_root_node->InvalidPointNum;
            if (father_ptr == STATIC_ROOT_NODE) rootNode = STATIC_ROOT_NODE->left_son_ptr;
            KDTreeNode * update_root = *Rebuild_Ptr;
            while (update_root != nullptr && update_root != rootNode){
                update_root = update_root->father_ptr;
                if (update_root->working_flag) break;
                if (update_root == update_root->father_ptr->left_son_ptr && update_root->father_ptr->need_push_down_to_left) break;
                if (update_root == update_root->father_ptr->right_son_ptr && update_root->father_ptr->need_push_down_to_right) break;
                Update(update_root);
            }
            pthread_mutex_lock(&search_flag_mutex);
            search_mutex_counter = 0;
            pthread_mutex_unlock(&search_flag_mutex);
            Rebuild_Ptr = nullptr;
            pthread_mutex_unlock(&working_flag_mutex);
            rebuildFlag = false;
            /* Delete discarded tree nodes */
            delete_tree_nodes(&old_root_node);
        } else {
            pthread_mutex_unlock(&working_flag_mutex);             
        }
        pthread_mutex_unlock(&rebuild_ptr_mutex_lock);         
        pthread_mutex_lock(&termination_flag_mutex_lock);
        terminated = terminationFlag;
        pthread_mutex_unlock(&termination_flag_mutex_lock);
        usleep(100); 
    }
    printf("rebuild thread terminated normally\n");    
}

template <typename PointT>
void KDTree<PointT>::run_operation(KDTreeNode ** root, OperationLoggerType operation){
    switch (operation.op)
    {
    case AddPoint:      
        point_add(root, operation.point, false, (*root)->division_axis);          
        break;
    case AddBox:
        range_add(root, operation.boxpoint, false);
        break;
    case DeletePoint:
        point_delete(root, operation.point, false);
        break;
    case DeleteBox:
        range_delete(root, operation.boxpoint, false, false);
        break;
    case DownsampleDelete:
        range_delete(root, operation.boxpoint, false, true);
        break;
    case PushDown:
        (*root)->tree_downsample_deleted |= operation.tree_downsample_deleted;
        (*root)->point_downsample_deleted |= operation.tree_downsample_deleted;
        (*root)->tree_deleted = operation.tree_deleted || (*root)->tree_downsample_deleted;
        (*root)->point_deleted = (*root)->tree_deleted || (*root)->point_downsample_deleted;
        if (operation.tree_downsample_deleted) (*root)->down_del_num = (*root)->TreeSize;
        if (operation.tree_deleted) (*root)->InvalidPointNum = (*root)->TreeSize;
            else (*root)->InvalidPointNum = (*root)->down_del_num;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;     
        break;
    default:
        break;
    }
}

template <typename PointT>
void KDTree<PointT>::build(PointVector point_cloud){
    if (rootNode != nullptr){
        delete_tree_nodes(&rootNode);
    }
    if (point_cloud.size() == 0) return;
    STATIC_ROOT_NODE = new KDTreeNode;
    init_tree_node(STATIC_ROOT_NODE); 
    buildTree(&STATIC_ROOT_NODE->left_son_ptr, 0, point_cloud.size()-1, point_cloud);
    Update(STATIC_ROOT_NODE);
    STATIC_ROOT_NODE->TreeSize = 0;
    rootNode = STATIC_ROOT_NODE->left_son_ptr;
}

template <typename PointT>
void KDTree<PointT>::nearest_search(PointT point, int k_nearest, PointVector& Nearest_Points, std::vector<float> & Point_Distance, double max_dist){
    Heap q(2*k_nearest);
    q.clear();
    std::vector<float> ().swap(Point_Distance);
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
        search(rootNode, k_nearest, point, q, max_dist);
    } else {
        pthread_mutex_lock(&search_flag_mutex);
        while (search_mutex_counter == -1)
        {
            pthread_mutex_unlock(&search_flag_mutex);
            usleep(1);
            pthread_mutex_lock(&search_flag_mutex);
        }
        search_mutex_counter += 1;
        pthread_mutex_unlock(&search_flag_mutex);  
        search(rootNode, k_nearest, point, q, max_dist);
        pthread_mutex_lock(&search_flag_mutex);
        search_mutex_counter -= 1;
        pthread_mutex_unlock(&search_flag_mutex);      
    }
    int k_found = min(k_nearest,int(q.size()));
    PointVector ().swap(Nearest_Points);
    std::vector<float> ().swap(Point_Distance);
    for (int i=0;i < k_found;i++){
        Nearest_Points.insert(Nearest_Points.begin(), q.top().point);
        Point_Distance.insert(Point_Distance.begin(), q.top().dist);
        q.pop();
    }
    return;
}

template <typename PointT>
void KDTree<PointT>::box_search(const BoxPointType &Box_of_Point, PointVector &Storage)
{
    Storage.clear();
    range_search(rootNode, Box_of_Point, Storage);
}

template <typename PointT>
void KDTree<PointT>::radius_search(PointT point, const float radius, PointVector &Storage)
{
    Storage.clear();
    radius_search(rootNode, point, radius, Storage);
}

template <typename PointT>
int KDTree<PointT>::add_points(PointVector & PointToAdd, bool downsample_on){
    int NewPointSize = PointToAdd.size();
    int tree_size = size();
    BoxPointType Box_of_Point;
    PointT downsample_result, mid_point;
    bool downsample_switch = downsample_on && DoDownsample;
    float min_dist, tmp_dist;
    int tmp_counter = 0;
    for (int i=0; i<PointToAdd.size();i++){
        if (downsample_switch){
            Box_of_Point.vertex_min[0] = floor(PointToAdd[i].x/downsampleSize)*downsampleSize;
            Box_of_Point.vertex_max[0] = Box_of_Point.vertex_min[0]+downsampleSize;
            Box_of_Point.vertex_min[1] = floor(PointToAdd[i].y/downsampleSize)*downsampleSize;
            Box_of_Point.vertex_max[1] = Box_of_Point.vertex_min[1]+downsampleSize;
            Box_of_Point.vertex_min[2] = floor(PointToAdd[i].z/downsampleSize)*downsampleSize;
            Box_of_Point.vertex_max[2] = Box_of_Point.vertex_min[2]+downsampleSize;
            mid_point.x = Box_of_Point.vertex_min[0] + (Box_of_Point.vertex_max[0]-Box_of_Point.vertex_min[0])/2.0;
            mid_point.y = Box_of_Point.vertex_min[1] + (Box_of_Point.vertex_max[1]-Box_of_Point.vertex_min[1])/2.0;
            mid_point.z = Box_of_Point.vertex_min[2] + (Box_of_Point.vertex_max[2]-Box_of_Point.vertex_min[2])/2.0;
            PointVector ().swap(downsampleStorage);
            range_search(rootNode, Box_of_Point, downsampleStorage);
            min_dist = calc_dist(PointToAdd[i],mid_point);
            downsample_result = PointToAdd[i];                
            for (int index = 0; index < downsampleStorage.size(); index++){
                tmp_dist = calc_dist(downsampleStorage[index], mid_point);
                if (tmp_dist < min_dist){
                    min_dist = tmp_dist;
                    downsample_result = downsampleStorage[index];
                }
            }
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
                if (downsampleStorage.size() > 1 || is_same_point(PointToAdd[i], downsample_result)){
                    if (downsampleStorage.size() > 0) range_delete(&rootNode, Box_of_Point, true, true);
                    point_add(&rootNode, downsample_result, true, rootNode->division_axis);
                    tmp_counter ++;                      
                }
            } else {
                if (downsampleStorage.size() > 1 || is_same_point(PointToAdd[i], downsample_result)){
                    OperationLoggerType  operation_delete, operation;
                    operation_delete.boxpoint = Box_of_Point;
                    operation_delete.op = DownsampleDelete;
                    operation.point = downsample_result;
                    operation.op = AddPoint;
                    pthread_mutex_lock(&working_flag_mutex);
                    if (downsampleStorage.size() > 0) range_delete(&rootNode, Box_of_Point, false , true);
                    point_add(&rootNode, downsample_result, false, rootNode->division_axis);
                    tmp_counter ++;
                    if (rebuildFlag){
                        pthread_mutex_lock(&rebuild_logger_mutex_lock);
                        if (downsampleStorage.size() > 0) rebuildLogger.push(operation_delete);
                        rebuildLogger.push(operation);
                        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                    }
                    pthread_mutex_unlock(&working_flag_mutex);
                };
            }
        } else {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
                point_add(&rootNode, PointToAdd[i], true, rootNode->division_axis);
            } else {
                OperationLoggerType operation;
                operation.point = PointToAdd[i];
                operation.op = AddPoint;                
                pthread_mutex_lock(&working_flag_mutex);
                point_add(&rootNode, PointToAdd[i], false, rootNode->division_axis);
                if (rebuildFlag){
                    pthread_mutex_lock(&rebuild_logger_mutex_lock);
                    rebuildLogger.push(operation);
                    pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                }
                pthread_mutex_unlock(&working_flag_mutex);       
            }
        }
    }
    return tmp_counter;
}

template <typename PointT>
void KDTree<PointT>::add_point_boxes(std::vector<BoxPointType> & BoxPoints){
    for (int i=0;i < BoxPoints.size();i++){
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
            range_add(&rootNode ,BoxPoints[i], true);
        } else {
            OperationLoggerType operation;
            operation.boxpoint = BoxPoints[i];
            operation.op = AddBox;
            pthread_mutex_lock(&working_flag_mutex);
            range_add(&rootNode ,BoxPoints[i], false);
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }               
            pthread_mutex_unlock(&working_flag_mutex);
        }    
    } 
    return;
}

template <typename PointT>
void KDTree<PointT>::delete_points(PointVector & PointToDel){
    for (int i=0;i<PointToDel.size();i++){
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
            point_delete(&rootNode, PointToDel[i], true);
        } else {
            OperationLoggerType operation;
            operation.point = PointToDel[i];
            operation.op = DeletePoint;
            pthread_mutex_lock(&working_flag_mutex);        
            point_delete(&rootNode, PointToDel[i], false);
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }      
    }      
    return;
}

template <typename PointT>
int KDTree<PointT>::delete_point_boxes(std::vector<BoxPointType> & BoxPoints){
    int tmp_counter = 0;
    for (int i=0;i < BoxPoints.size();i++){ 
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != rootNode){
            tmp_counter += range_delete(&rootNode ,BoxPoints[i], true, false);
        } else {
            OperationLoggerType operation;
            operation.boxpoint = BoxPoints[i];
            operation.op = DeleteBox;     
            pthread_mutex_lock(&working_flag_mutex); 
            tmp_counter += range_delete(&rootNode ,BoxPoints[i], false, false);
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }                
            pthread_mutex_unlock(&working_flag_mutex);
        }
    } 
    return tmp_counter;
}

template <typename PointT>
void KDTree<PointT>::acquire_removed_points(PointVector & removed_points){
    pthread_mutex_lock(&points_deleted_rebuild_mutex_lock); 
    for (int i = 0; i < deletedPoints.size();i++){
        removed_points.push_back(deletedPoints[i]);
    }
    for (int i = 0; i < deletedPointsMultithread.size();i++){
        removed_points.push_back(deletedPointsMultithread[i]);
    }
    deletedPoints.clear();
    deletedPointsMultithread.clear();
    pthread_mutex_unlock(&points_deleted_rebuild_mutex_lock);   
    return;
}

template <typename PointT>
void KDTree<PointT>::buildTree(KDTreeNode ** root, int l, int r, PointVector & Storage){
    if (l>r) return;
    *root = new KDTreeNode;
    init_tree_node(*root);
    int mid = (l+r)>>1;
    int div_axis = 0;
    int i;
    // Find the best division Axis
    float min_value[3] = {INFINITY, INFINITY, INFINITY};
    float max_value[3] = {-INFINITY, -INFINITY, -INFINITY};
    float dim_range[3] = {0,0,0};
    for (i=l;i<=r;i++){
        min_value[0] = min(min_value[0], Storage[i].x);
        min_value[1] = min(min_value[1], Storage[i].y);
        min_value[2] = min(min_value[2], Storage[i].z);
        max_value[0] = max(max_value[0], Storage[i].x);
        max_value[1] = max(max_value[1], Storage[i].y);
        max_value[2] = max(max_value[2], Storage[i].z);
    }
    // Select the longest dimension as division axis
    for (i=0;i<3;i++) dim_range[i] = max_value[i] - min_value[i];
    for (i=1;i<3;i++) if (dim_range[i] > dim_range[div_axis]) div_axis = i;
    // Divide by the division axis and recursively build.

    (*root)->division_axis = div_axis;
    switch (div_axis)
    {
    case 0:
        nth_element(begin(Storage)+l, begin(Storage)+mid, begin(Storage)+r+1, point_cmp_x);
        break;
    case 1:
        nth_element(begin(Storage)+l, begin(Storage)+mid, begin(Storage)+r+1, point_cmp_y);
        break;
    case 2:
        nth_element(begin(Storage)+l, begin(Storage)+mid, begin(Storage)+r+1, point_cmp_z);
        break;
    default:
        nth_element(begin(Storage)+l, begin(Storage)+mid, begin(Storage)+r+1, point_cmp_x);
        break;
    }  
    (*root)->point = Storage[mid]; 
    KDTreeNode * left_son = nullptr, * right_son = nullptr;
    buildTree(&left_son, l, mid-1, Storage);
    buildTree(&right_son, mid+1, r, Storage);  
    (*root)->left_son_ptr = left_son;
    (*root)->right_son_ptr = right_son;
    Update((*root));  
    return;
}

template <typename PointT>
void KDTree<PointT>::rebuild(KDTreeNode ** root){    
    KDTreeNode * father_ptr;
    if ((*root)->TreeSize >= MultiThreadRebuildPointNum) { 
        if (!pthread_mutex_trylock(&rebuild_ptr_mutex_lock)){     
            if (Rebuild_Ptr == nullptr || ((*root)->TreeSize > (*Rebuild_Ptr)->TreeSize)) {
                Rebuild_Ptr = root;          
            }
            pthread_mutex_unlock(&rebuild_ptr_mutex_lock);
        }
    } else {
        father_ptr = (*root)->father_ptr;
        int size_rec = (*root)->TreeSize;
        PCLStorage.clear();
        flatten(*root, PCLStorage, DeletePointsRecord);
        delete_tree_nodes(root);
        buildTree(root, 0, PCLStorage.size()-1, PCLStorage);
        if (*root != nullptr) (*root)->father_ptr = father_ptr;
        if (*root == rootNode) STATIC_ROOT_NODE->left_son_ptr = *root;
    } 
    return;
}

template <typename PointT>
int KDTree<PointT>::range_delete(KDTreeNode ** root,  BoxPointType boxpoint, bool allow_rebuild, bool is_downsample){   
    if ((*root) == nullptr || (*root)->tree_deleted) return 0;
    (*root)->working_flag = true;
    push_down(*root);
    int tmp_counter = 0;
    if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] || boxpoint.vertex_min[0] > (*root)->node_range_x[1]) return 0;
    if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] || boxpoint.vertex_min[1] > (*root)->node_range_y[1]) return 0;
    if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] || boxpoint.vertex_min[2] > (*root)->node_range_z[1]) return 0;
    if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] && boxpoint.vertex_max[0] > (*root)->node_range_x[1] && boxpoint.vertex_min[1] <= (*root)->node_range_y[0] && boxpoint.vertex_max[1] > (*root)->node_range_y[1] && boxpoint.vertex_min[2] <= (*root)->node_range_z[0] && boxpoint.vertex_max[2] > (*root)->node_range_z[1]){
        (*root)->tree_deleted = true;
        (*root)->point_deleted = true;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;
        tmp_counter = (*root)->TreeSize - (*root)->InvalidPointNum;
        (*root)->InvalidPointNum = (*root)->TreeSize;
        if (is_downsample){
            (*root)->tree_downsample_deleted = true;
            (*root)->point_downsample_deleted = true;
            (*root)->down_del_num = (*root)->TreeSize;
        }
        return tmp_counter;
    }
    if (!(*root)->point_deleted && boxpoint.vertex_min[0] <= (*root)->point.x && boxpoint.vertex_max[0] > (*root)->point.x && boxpoint.vertex_min[1] <= (*root)->point.y && boxpoint.vertex_max[1] > (*root)->point.y && boxpoint.vertex_min[2] <= (*root)->point.z && boxpoint.vertex_max[2] > (*root)->point.z){
        (*root)->point_deleted = true;
        tmp_counter += 1;
        if (is_downsample) (*root)->point_downsample_deleted = true;
    }
    OperationLoggerType delete_box_log;
    struct timespec Timeout;    
    if (is_downsample) delete_box_log.op = DownsampleDelete;
        else delete_box_log.op = DeleteBox;
    delete_box_log.boxpoint = boxpoint;
    if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr){
        tmp_counter += range_delete(&((*root)->left_son_ptr), boxpoint, allow_rebuild, is_downsample);
    } else {
        pthread_mutex_lock(&working_flag_mutex);
        tmp_counter += range_delete(&((*root)->left_son_ptr), boxpoint, false, is_downsample);
        if (rebuildFlag){
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            rebuildLogger.push(delete_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }
    if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr){
        tmp_counter += range_delete(&((*root)->right_son_ptr), boxpoint, allow_rebuild, is_downsample);
    } else {
        pthread_mutex_lock(&working_flag_mutex);
        tmp_counter += range_delete(&((*root)->right_son_ptr), boxpoint, false, is_downsample);
        if (rebuildFlag){
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            rebuildLogger.push(delete_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }    
    Update(*root);     
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < MultiThreadRebuildPointNum) Rebuild_Ptr = nullptr; 
    bool need_rebuild = allow_rebuild & criterion_check((*root));
    if (need_rebuild) rebuild(root);
    if ((*root) != nullptr) (*root)->working_flag = false;
    return tmp_counter;
}

template <typename PointT>
void KDTree<PointT>::point_delete(KDTreeNode ** root, PointT point, bool allow_rebuild){   
    if ((*root) == nullptr || (*root)->tree_deleted) return;
    (*root)->working_flag = true;
    push_down(*root);
    if (is_same_point((*root)->point, point) && !(*root)->point_deleted) {
        (*root)->point_deleted = true;
        (*root)->InvalidPointNum += 1;
        if ((*root)->InvalidPointNum == (*root)->TreeSize) (*root)->tree_deleted = true;    
        return;
    }
    OperationLoggerType delete_log;
    struct timespec Timeout;    
    delete_log.op = DeletePoint;
    delete_log.point = point;     
    if (((*root)->division_axis == 0 && point.x < (*root)->point.x) || ((*root)->division_axis == 1 && point.y < (*root)->point.y) || ((*root)->division_axis == 2 && point.z < (*root)->point.z)){           
        if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr){          
            point_delete(&(*root)->left_son_ptr, point, allow_rebuild);         
        } else {
            pthread_mutex_lock(&working_flag_mutex);
            point_delete(&(*root)->left_son_ptr, point,false);
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(delete_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    } else {       
        if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr){         
            point_delete(&(*root)->right_son_ptr, point, allow_rebuild);         
        } else {
            pthread_mutex_lock(&working_flag_mutex); 
            point_delete(&(*root)->right_son_ptr, point, false);
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(delete_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }        
    }
    Update(*root);
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < MultiThreadRebuildPointNum) Rebuild_Ptr = nullptr; 
    bool need_rebuild = allow_rebuild & criterion_check((*root));
    if (need_rebuild) rebuild(root);
    if ((*root) != nullptr) (*root)->working_flag = false;   
    return;
}

template <typename PointT>
void KDTree<PointT>::range_add(KDTreeNode ** root, BoxPointType boxpoint, bool allow_rebuild){
    if ((*root) == nullptr) return;
    (*root)->working_flag = true;
    push_down(*root);
    if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] || boxpoint.vertex_min[0] > (*root)->node_range_x[1]) return;
    if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] || boxpoint.vertex_min[1] > (*root)->node_range_y[1]) return;
    if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] || boxpoint.vertex_min[2] > (*root)->node_range_z[1]) return;
    if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] && boxpoint.vertex_max[0] > (*root)->node_range_x[1] && boxpoint.vertex_min[1] <= (*root)->node_range_y[0] && boxpoint.vertex_max[1]> (*root)->node_range_y[1] && boxpoint.vertex_min[2] <= (*root)->node_range_z[0] && boxpoint.vertex_max[2] > (*root)->node_range_z[1]){
        (*root)->tree_deleted = false || (*root)->tree_downsample_deleted;
        (*root)->point_deleted = false || (*root)->point_downsample_deleted;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;
        (*root)->InvalidPointNum = (*root)->down_del_num; 
        return;
    }
    if (boxpoint.vertex_min[0] <= (*root)->point.x && boxpoint.vertex_max[0] > (*root)->point.x && boxpoint.vertex_min[1] <= (*root)->point.y && boxpoint.vertex_max[1] > (*root)->point.y && boxpoint.vertex_min[2] <= (*root)->point.z && boxpoint.vertex_max[2] > (*root)->point.z){
        (*root)->point_deleted = (*root)->point_downsample_deleted;
    }
    OperationLoggerType add_box_log;
    struct timespec Timeout;    
    add_box_log.op = AddBox;
    add_box_log.boxpoint = boxpoint;
    if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr){
        range_add(&((*root)->left_son_ptr), boxpoint, allow_rebuild);
    } else {
        pthread_mutex_lock(&working_flag_mutex);
        range_add(&((*root)->left_son_ptr), boxpoint, false);
        if (rebuildFlag){
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            rebuildLogger.push(add_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
        }        
        pthread_mutex_unlock(&working_flag_mutex);
    }
    if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr){
        range_add(&((*root)->right_son_ptr), boxpoint, allow_rebuild);
    } else {
        pthread_mutex_lock(&working_flag_mutex);
        range_add(&((*root)->right_son_ptr), boxpoint, false);
        if (rebuildFlag){
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            rebuildLogger.push(add_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }
    Update(*root);
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < MultiThreadRebuildPointNum) Rebuild_Ptr = nullptr; 
    bool need_rebuild = allow_rebuild & criterion_check((*root));
    if (need_rebuild) rebuild(root);
    if ((*root) != nullptr) (*root)->working_flag = false;   
    return;
}

template <typename PointT>
void KDTree<PointT>::point_add(KDTreeNode ** root, PointT point, bool allow_rebuild, int father_axis){     
    if (*root == nullptr){
        *root = new KDTreeNode;
        init_tree_node(*root);
        (*root)->point = point;
        (*root)->division_axis = (father_axis + 1) % 3;
        Update(*root);
        return;
    }
    (*root)->working_flag = true;
    OperationLoggerType add_log;
    struct timespec Timeout;    
    add_log.op = AddPoint;
    add_log.point = point;
    push_down(*root);
    if (((*root)->division_axis == 0 && point.x < (*root)->point.x) || ((*root)->division_axis == 1 && point.y < (*root)->point.y) || ((*root)->division_axis == 2 && point.z < (*root)->point.z)){
        if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr){          
            point_add(&(*root)->left_son_ptr, point, allow_rebuild, (*root)->division_axis);
        } else {
            pthread_mutex_lock(&working_flag_mutex);
            point_add(&(*root)->left_son_ptr, point, false,(*root)->division_axis);
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(add_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
            }
            pthread_mutex_unlock(&working_flag_mutex);            
        }
    } else {  
        if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr){         
            point_add(&(*root)->right_son_ptr, point, allow_rebuild,(*root)->division_axis);
        } else {
            pthread_mutex_lock(&working_flag_mutex);
            point_add(&(*root)->right_son_ptr, point, false,(*root)->division_axis);       
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(add_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);                 
            }
            pthread_mutex_unlock(&working_flag_mutex); 
        }
    }
    Update(*root);   
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < MultiThreadRebuildPointNum) Rebuild_Ptr = nullptr; 
    bool need_rebuild = allow_rebuild & criterion_check((*root));
    if (need_rebuild) rebuild(root); 
    if ((*root) != nullptr) (*root)->working_flag = false;   
    return;
}

template <typename PointT>
void KDTree<PointT>::search(KDTreeNode * root, int k_nearest, PointT point, Heap &q, double max_dist){
    if (root == nullptr || root->tree_deleted) return;   
    double cur_dist = calc_box_dist(root, point);
    double max_dist_sqr = max_dist * max_dist;
    if (cur_dist > max_dist_sqr) return;    
    int retval; 
    if (root->need_push_down_to_left || root->need_push_down_to_right) {
        retval = pthread_mutex_trylock(&(root->push_down_mutex_lock));
        if (retval == 0){
            push_down(root);
            pthread_mutex_unlock(&(root->push_down_mutex_lock));
        } else {
            pthread_mutex_lock(&(root->push_down_mutex_lock));
            pthread_mutex_unlock(&(root->push_down_mutex_lock));
        }
    }
    if (!root->point_deleted){
        float dist = calc_dist(point, root->point);
        if (dist <= max_dist_sqr && (q.size() < k_nearest || dist < q.top().dist)){
            if (q.size() >= k_nearest) q.pop();
            HeapPoint current_point{root->point, dist};                    
            q.push(current_point);            
        }
    }  
    int cur_search_counter;
    float dist_left_node = calc_box_dist(root->left_son_ptr, point);
    float dist_right_node = calc_box_dist(root->right_son_ptr, point);
    if (q.size()< k_nearest || dist_left_node < q.top().dist && dist_right_node < q.top().dist){
        if (dist_left_node <= dist_right_node) {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr){
                search(root->left_son_ptr, k_nearest, point, q, max_dist);
            } else {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1)
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);
                search(root->left_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
            if (q.size() < k_nearest || dist_right_node < q.top().dist) {
                if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr){
                    search(root->right_son_ptr, k_nearest, point, q, max_dist);
                } else {
                    pthread_mutex_lock(&search_flag_mutex);
                    while (search_mutex_counter == -1)
                    {
                        pthread_mutex_unlock(&search_flag_mutex);
                        usleep(1);
                        pthread_mutex_lock(&search_flag_mutex);
                    }
                    search_mutex_counter += 1;
                    pthread_mutex_unlock(&search_flag_mutex);                    
                    search(root->right_son_ptr, k_nearest, point, q, max_dist);
                    pthread_mutex_lock(&search_flag_mutex);
                    search_mutex_counter -= 1;
                    pthread_mutex_unlock(&search_flag_mutex);
                }                
            }
        } else {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr){
                search(root->right_son_ptr, k_nearest, point, q, max_dist);
            } else {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1)
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);                   
                search(root->right_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
            if (q.size() < k_nearest || dist_left_node < q.top().dist) {            
                if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr){
                    search(root->left_son_ptr, k_nearest, point, q, max_dist);
                } else {
                    pthread_mutex_lock(&search_flag_mutex);
                    while (search_mutex_counter == -1)
                    {
                        pthread_mutex_unlock(&search_flag_mutex);
                        usleep(1);
                        pthread_mutex_lock(&search_flag_mutex);
                    }
                    search_mutex_counter += 1;
                    pthread_mutex_unlock(&search_flag_mutex);  
                    search(root->left_son_ptr, k_nearest, point, q, max_dist);
                    pthread_mutex_lock(&search_flag_mutex);
                    search_mutex_counter -= 1;
                    pthread_mutex_unlock(&search_flag_mutex);
                }
            }
        }
    } else {
        if (dist_left_node < q.top().dist) {        
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr){
                search(root->left_son_ptr, k_nearest, point, q, max_dist);
            } else {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1)
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);  
                search(root->left_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
        }
        if (dist_right_node < q.top().dist) {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr){
                search(root->right_son_ptr, k_nearest, point, q, max_dist);
            } else {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1)
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);  
                search(root->right_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
        }
    }
    return;
}

template <typename PointT>
void KDTree<PointT>::range_search(KDTreeNode *root, BoxPointType boxpoint, PointVector & Storage){
    if (root == nullptr) return;
    push_down(root);
    if (boxpoint.vertex_max[0] <= root->node_range_x[0] || boxpoint.vertex_min[0] > root->node_range_x[1]) return;
    if (boxpoint.vertex_max[1] <= root->node_range_y[0] || boxpoint.vertex_min[1] > root->node_range_y[1]) return;
    if (boxpoint.vertex_max[2] <= root->node_range_z[0] || boxpoint.vertex_min[2] > root->node_range_z[1]) return;
    if (boxpoint.vertex_min[0] <= root->node_range_x[0] && boxpoint.vertex_max[0] > root->node_range_x[1] && boxpoint.vertex_min[1] <= root->node_range_y[0] && boxpoint.vertex_max[1] > root->node_range_y[1] && boxpoint.vertex_min[2] <= root->node_range_z[0] && boxpoint.vertex_max[2] > root->node_range_z[1]){
        flatten(root, Storage, NotRecord);
        return;
    }
    if (boxpoint.vertex_min[0] <= root->point.x && boxpoint.vertex_max[0] > root->point.x && boxpoint.vertex_min[1] <= root->point.y && boxpoint.vertex_max[1] > root->point.y && boxpoint.vertex_min[2] <= root->point.z && boxpoint.vertex_max[2] > root->point.z){
        if (!root->point_deleted) Storage.push_back(root->point);
    }
    if ((Rebuild_Ptr == nullptr) || root->left_son_ptr != *Rebuild_Ptr){
        range_search(root->left_son_ptr, boxpoint, Storage);
    } else {
        pthread_mutex_lock(&search_flag_mutex);
        range_search(root->left_son_ptr, boxpoint, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }
    if ((Rebuild_Ptr == nullptr) || root->right_son_ptr != *Rebuild_Ptr){
        range_search(root->right_son_ptr, boxpoint, Storage);
    } else {
        pthread_mutex_lock(&search_flag_mutex);
        range_search(root->right_son_ptr, boxpoint, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }
    return;    
}

template <typename PointT>
void KDTree<PointT>::radius_search(KDTreeNode *root, PointT point, float radius, PointVector &Storage)
{
    if (root == nullptr)
        return;
    push_down(root);
    PointT range_center;
    range_center.x = (root->node_range_x[0] + root->node_range_x[1]) * 0.5;
    range_center.y = (root->node_range_y[0] + root->node_range_y[1]) * 0.5;
    range_center.z = (root->node_range_z[0] + root->node_range_z[1]) * 0.5;
    float dist = sqrt(calc_dist(range_center, point));
    if (dist > radius + sqrt(root->radius_sq)) return;
    if (dist <= radius - sqrt(root->radius_sq)) 
    {
        flatten(root, Storage, NotRecord);
        return;
    }
    if (!root->point_deleted && calc_dist(root->point, point) <= radius * radius){
        Storage.push_back(root->point);
    }
    if ((Rebuild_Ptr == nullptr) || root->left_son_ptr != *Rebuild_Ptr)
    {
        radius_search(root->left_son_ptr, point, radius, Storage);
    }
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        radius_search(root->left_son_ptr, point, radius, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }
    if ((Rebuild_Ptr == nullptr) || root->right_son_ptr != *Rebuild_Ptr)
    {
        radius_search(root->right_son_ptr, point, radius, Storage);
    }
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        radius_search(root->right_son_ptr, point, radius, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }    
    return;
}

template <typename PointT>
bool KDTree<PointT>::criterion_check(KDTreeNode * root){
    if (root->TreeSize <= MinUnbalancedTreeSize){
        return false;
    }
    float balance_evaluation = 0.0f;
    float delete_evaluation = 0.0f;
    KDTreeNode * son_ptr = root->left_son_ptr;
    if (son_ptr == nullptr) son_ptr = root->right_son_ptr;
    delete_evaluation = float(root->InvalidPointNum)/ root->TreeSize;
    balance_evaluation = float(son_ptr->TreeSize) / (root->TreeSize-1);  
    if (delete_evaluation > deleteParam){
        return true;
    }
    if (balance_evaluation > balanceParam || balance_evaluation < 1-balanceParam){
        return true;
    } 
    return false;
}

template <typename PointT>
void KDTree<PointT>::push_down(KDTreeNode *root){
    if (root == nullptr) return;
    OperationLoggerType operation;
    operation.op = PushDown;
    operation.tree_deleted = root->tree_deleted;
    operation.tree_downsample_deleted = root->tree_downsample_deleted;
    if (root->need_push_down_to_left && root->left_son_ptr != nullptr){
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr){
            root->left_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->tree_deleted = root->tree_deleted || root->left_son_ptr->tree_downsample_deleted;
            root->left_son_ptr->point_deleted = root->left_son_ptr->tree_deleted || root->left_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted) root->left_son_ptr->down_del_num = root->left_son_ptr->TreeSize;
            if (root->tree_deleted) root->left_son_ptr->InvalidPointNum = root->left_son_ptr->TreeSize;
                else root->left_son_ptr->InvalidPointNum = root->left_son_ptr->down_del_num;
            root->left_son_ptr->need_push_down_to_left = true;
            root->left_son_ptr->need_push_down_to_right = true;
            root->need_push_down_to_left = false;                
        } else {
            pthread_mutex_lock(&working_flag_mutex);
            root->left_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->tree_deleted = root->tree_deleted || root->left_son_ptr->tree_downsample_deleted;
            root->left_son_ptr->point_deleted = root->left_son_ptr->tree_deleted || root->left_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted) root->left_son_ptr->down_del_num = root->left_son_ptr->TreeSize;
            if (root->tree_deleted) root->left_son_ptr->InvalidPointNum = root->left_son_ptr->TreeSize;
                else root->left_son_ptr->InvalidPointNum = root->left_son_ptr->down_del_num;            
            root->left_son_ptr->need_push_down_to_left = true;
            root->left_son_ptr->need_push_down_to_right = true;
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            root->need_push_down_to_left = false;
            pthread_mutex_unlock(&working_flag_mutex);            
        }
    }
    if (root->need_push_down_to_right && root->right_son_ptr != nullptr){
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr){
            root->right_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->tree_deleted = root->tree_deleted || root->right_son_ptr->tree_downsample_deleted;
            root->right_son_ptr->point_deleted = root->right_son_ptr->tree_deleted || root->right_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted) root->right_son_ptr->down_del_num = root->right_son_ptr->TreeSize;
            if (root->tree_deleted) root->right_son_ptr->InvalidPointNum = root->right_son_ptr->TreeSize;
                else root->right_son_ptr->InvalidPointNum = root->right_son_ptr->down_del_num;
            root->right_son_ptr->need_push_down_to_left = true;
            root->right_son_ptr->need_push_down_to_right = true;
            root->need_push_down_to_right = false;
        } else {
            pthread_mutex_lock(&working_flag_mutex);
            root->right_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->tree_deleted = root->tree_deleted || root->right_son_ptr->tree_downsample_deleted;
            root->right_son_ptr->point_deleted = root->right_son_ptr->tree_deleted || root->right_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted) root->right_son_ptr->down_del_num = root->right_son_ptr->TreeSize;
            if (root->tree_deleted) root->right_son_ptr->InvalidPointNum = root->right_son_ptr->TreeSize;
                else root->right_son_ptr->InvalidPointNum = root->right_son_ptr->down_del_num;            
            root->right_son_ptr->need_push_down_to_left = true;
            root->right_son_ptr->need_push_down_to_right = true;
            if (rebuildFlag){
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                rebuildLogger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }            
            root->need_push_down_to_right = false;
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    return;
}

template <typename PointT>
void KDTree<PointT>::Update(KDTreeNode * root){
    KDTreeNode * left_son_ptr = root->left_son_ptr;
    KDTreeNode * right_son_ptr = root->right_son_ptr;
    float tmp_range_x[2] = {INFINITY, -INFINITY};
    float tmp_range_y[2] = {INFINITY, -INFINITY};
    float tmp_range_z[2] = {INFINITY, -INFINITY};
    // Update Tree Size   
    if (left_son_ptr != nullptr && right_son_ptr != nullptr){
        root->TreeSize = left_son_ptr->TreeSize + right_son_ptr->TreeSize + 1;
        root->InvalidPointNum = left_son_ptr->InvalidPointNum + right_son_ptr->InvalidPointNum + (root->point_deleted? 1:0);
        root->down_del_num = left_son_ptr->down_del_num + right_son_ptr->down_del_num + (root->point_downsample_deleted? 1:0);
        root->tree_downsample_deleted = left_son_ptr->tree_downsample_deleted & right_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
        root->tree_deleted = left_son_ptr->tree_deleted && right_son_ptr->tree_deleted && root->point_deleted;
        if (root->tree_deleted || (!left_son_ptr->tree_deleted && !right_son_ptr->tree_deleted && !root->point_deleted)){
            tmp_range_x[0] = min(min(left_son_ptr->node_range_x[0],right_son_ptr->node_range_x[0]),root->point.x);
            tmp_range_x[1] = max(max(left_son_ptr->node_range_x[1],right_son_ptr->node_range_x[1]),root->point.x);
            tmp_range_y[0] = min(min(left_son_ptr->node_range_y[0],right_son_ptr->node_range_y[0]),root->point.y);
            tmp_range_y[1] = max(max(left_son_ptr->node_range_y[1],right_son_ptr->node_range_y[1]),root->point.y);
            tmp_range_z[0] = min(min(left_son_ptr->node_range_z[0],right_son_ptr->node_range_z[0]),root->point.z);
            tmp_range_z[1] = max(max(left_son_ptr->node_range_z[1],right_son_ptr->node_range_z[1]),root->point.z);
        } else {
            if (!left_son_ptr->tree_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
            }
            if (!right_son_ptr->tree_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], right_son_ptr->node_range_z[1]);                
            }
            if (!root->point_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
                tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
                tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
                tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
                tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
                tmp_range_z[1] = max(tmp_range_z[1], root->point.z);                 
            }
        }
    } else if (left_son_ptr != nullptr){
        root->TreeSize = left_son_ptr->TreeSize + 1;
        root->InvalidPointNum = left_son_ptr->InvalidPointNum + (root->point_deleted?1:0);
        root->down_del_num = left_son_ptr->down_del_num + (root->point_downsample_deleted?1:0);
        root->tree_downsample_deleted = left_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
        root->tree_deleted = left_son_ptr->tree_deleted && root->point_deleted;
        if (root->tree_deleted || (!left_son_ptr->tree_deleted && !root->point_deleted)){
            tmp_range_x[0] = min(left_son_ptr->node_range_x[0],root->point.x);
            tmp_range_x[1] = max(left_son_ptr->node_range_x[1],root->point.x);
            tmp_range_y[0] = min(left_son_ptr->node_range_y[0],root->point.y);
            tmp_range_y[1] = max(left_son_ptr->node_range_y[1],root->point.y); 
            tmp_range_z[0] = min(left_son_ptr->node_range_z[0],root->point.z);
            tmp_range_z[1] = max(left_son_ptr->node_range_z[1],root->point.z);  
        } else {
            if (!left_son_ptr->tree_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], left_son_ptr->node_range_z[1]);                
            }
            if (!root->point_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
                tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
                tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
                tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
                tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
                tmp_range_z[1] = max(tmp_range_z[1], root->point.z);                 
            }            
        }

    } else if (right_son_ptr != nullptr){
        root->TreeSize = right_son_ptr->TreeSize + 1;
        root->InvalidPointNum = right_son_ptr->InvalidPointNum + (root->point_deleted? 1:0);
        root->down_del_num = right_son_ptr->down_del_num + (root->point_downsample_deleted? 1:0);        
        root->tree_downsample_deleted = right_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
        root->tree_deleted = right_son_ptr->tree_deleted && root->point_deleted;
        if (root->tree_deleted || (!right_son_ptr->tree_deleted && !root->point_deleted)){
            tmp_range_x[0] = min(right_son_ptr->node_range_x[0],root->point.x);
            tmp_range_x[1] = max(right_son_ptr->node_range_x[1],root->point.x);
            tmp_range_y[0] = min(right_son_ptr->node_range_y[0],root->point.y);
            tmp_range_y[1] = max(right_son_ptr->node_range_y[1],root->point.y); 
            tmp_range_z[0] = min(right_son_ptr->node_range_z[0],root->point.z);
            tmp_range_z[1] = max(right_son_ptr->node_range_z[1],root->point.z); 
        } else {
            if (!right_son_ptr->tree_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], right_son_ptr->node_range_z[1]);                
            }
            if (!root->point_deleted){
                tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
                tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
                tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
                tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
                tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
                tmp_range_z[1] = max(tmp_range_z[1], root->point.z);                 
            }            
        }
    } else {
        root->TreeSize = 1;
        root->InvalidPointNum = (root->point_deleted? 1:0);
        root->down_del_num = (root->point_downsample_deleted? 1:0);
        root->tree_downsample_deleted = root->point_downsample_deleted;
        root->tree_deleted = root->point_deleted;
        tmp_range_x[0] = root->point.x;
        tmp_range_x[1] = root->point.x;        
        tmp_range_y[0] = root->point.y;
        tmp_range_y[1] = root->point.y; 
        tmp_range_z[0] = root->point.z;
        tmp_range_z[1] = root->point.z;                 
    }
    memcpy(root->node_range_x,tmp_range_x,sizeof(tmp_range_x));
    memcpy(root->node_range_y,tmp_range_y,sizeof(tmp_range_y));
    memcpy(root->node_range_z,tmp_range_z,sizeof(tmp_range_z));
    float x_L = (root->node_range_x[1] - root->node_range_x[0]) * 0.5;
    float y_L = (root->node_range_y[1] - root->node_range_y[0]) * 0.5;
    float z_L = (root->node_range_z[1] - root->node_range_z[0]) * 0.5;
    root->radius_sq = x_L*x_L + y_L * y_L + z_L * z_L;    
    if (left_son_ptr != nullptr) left_son_ptr -> father_ptr = root;
    if (right_son_ptr != nullptr) right_son_ptr -> father_ptr = root;
    if (root == rootNode && root->TreeSize > 3){
        KDTreeNode * son_ptr = root->left_son_ptr;
        if (son_ptr == nullptr) son_ptr = root->right_son_ptr;
        float tmp_bal = float(son_ptr->TreeSize) / (root->TreeSize-1);
        root->alpha_del = float(root->InvalidPointNum)/ root->TreeSize;
        root->alpha_bal = (tmp_bal>=0.5-EPSS)?tmp_bal:1-tmp_bal;
    }   
    return;
}

template <typename PointT>
void KDTree<PointT>::flatten(KDTreeNode * root, PointVector &Storage, DeletePointStorageOpts storage_type){
    if (root == nullptr) return;
    push_down(root);
    if (!root->point_deleted) {
        Storage.push_back(root->point);
    }
    flatten(root->left_son_ptr, Storage, storage_type);
    flatten(root->right_son_ptr, Storage, storage_type);
    switch (storage_type)
    {
    case NotRecord:
        break;
    case DeletePointsRecord:
        if (root->point_deleted && !root->point_downsample_deleted) {
            deletedPoints.push_back(root->point);
        }       
        break;
    case MultiThreadRecord:
        if (root->point_deleted  && !root->point_downsample_deleted) {
            deletedPointsMultithread.push_back(root->point);
        }
        break;
    default:
        break;
    }     
    return;
}

template <typename PointT>
void KDTree<PointT>::delete_tree_nodes(KDTreeNode ** root){ 
    if (*root == nullptr) return;
    push_down(*root);
    delete_tree_nodes(&(*root)->left_son_ptr);
    delete_tree_nodes(&(*root)->right_son_ptr);
    
    pthread_mutex_destroy( &(*root)->push_down_mutex_lock);         
    delete *root;
    *root = nullptr;                    
}

template <typename PointT>
bool KDTree<PointT>::is_same_point(PointT a, PointT b){
    return (fabs(a.x-b.x) < EPSS && fabs(a.y-b.y) < EPSS && fabs(a.z-b.z) < EPSS );
}

template <typename PointT>
float KDTree<PointT>::calc_dist(PointT a, PointT b){
    float dist = 0.0f;
    dist = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);
    return dist;
}

template <typename PointT>
float KDTree<PointT>::calc_box_dist(KDTreeNode * node, PointT point){
    if (node == nullptr) return INFINITY;
    float min_dist = 0.0;
    if (point.x < node->node_range_x[0]) min_dist += (point.x - node->node_range_x[0])*(point.x - node->node_range_x[0]);
    if (point.x > node->node_range_x[1]) min_dist += (point.x - node->node_range_x[1])*(point.x - node->node_range_x[1]);
    if (point.y < node->node_range_y[0]) min_dist += (point.y - node->node_range_y[0])*(point.y - node->node_range_y[0]);
    if (point.y > node->node_range_y[1]) min_dist += (point.y - node->node_range_y[1])*(point.y - node->node_range_y[1]);
    if (point.z < node->node_range_z[0]) min_dist += (point.z - node->node_range_z[0])*(point.z - node->node_range_z[0]);
    if (point.z > node->node_range_z[1]) min_dist += (point.z - node->node_range_z[1])*(point.z - node->node_range_z[1]);
    return min_dist;
}

template <typename PointT> bool KDTree<PointT>::point_cmp_x(PointT a, PointT b) { return a.x < b.x;}
template <typename PointT> bool KDTree<PointT>::point_cmp_y(PointT a, PointT b) { return a.y < b.y;}
template <typename PointT> bool KDTree<PointT>::point_cmp_z(PointT a, PointT b) { return a.z < b.z;}

// manual queue
template <typename T>
void Queue<T>::clear(){
    head = 0;
    tail = 0;
    counter = 0;
    is_empty = true;
}

template <typename T>
void Queue<T>::pop(){
    if (counter == 0) return;
    head ++;
    head %= QSize;
    counter --;
    if (counter == 0) is_empty = true;
}

template <typename T>
T Queue<T>::front(){
    return q[head];
}

template <typename T>
T Queue<T>::back(){
    return q[tail];
}

template <typename T>
void Queue<T>::push(T op){
    q[tail] = op;
    counter ++;
    if (is_empty) is_empty = false;
    tail ++;
    tail %= QSize;
}

template <typename T>
bool Queue<T>::empty(){
    return is_empty;
}

template <typename T>
int Queue<T>::size(){
    return counter;
}

template class KDTree<pcl::PointXYZ>;
template class KDTree<pcl::PointXYZI>;
template class KDTree<pcl::PointXYZINormal>;