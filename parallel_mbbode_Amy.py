#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created on 2016.11.20
@author: Amy
'''
import time
import random
import math
from pyspark import SparkContext

def range2rect(size, num_var):
    #构造size x num_var 的矩阵
    res = range(size)
    for i in range(size):
        res[i] = range(num_var)
    return res

def make_population(size, num_var, rp, rm, f, p_mutate, time_base, lambdaa):
    #构造一个population，包含size个候选解chrom,每个chrom是由num_var个记录各vm所放置的pm下标组成的
    population = {
        'size': size,                                                       # population中chrom的个数
        'num_var': num_var,                                                 # chrom元素个数，即vm个数
        'rp': rp,                                                           # num_var个vm分别对应的cpu demand
        'rm': rm,
        'f': f,                                                             # 差分因子，系统运行前必须指定
        'p_mutate': p_mutate,                                               # 突变概率
        'time_base': time_base,                                             # 时间基数
        'lambdaa': tuple(lambdaa),                                          # 每个chrom被选中作为迁入岛屿的概率，rank值越大，概率越大
        'population': range2rect(size, num_var),                            # size x num_var的矩阵
        'population0': range2rect(size, num_var),                           # 用以保存初始化的chroms
        'p_cost': range2rect(size, num_var),                                # 记录num_var个pm被使用的cpu量 size x num_var的矩阵
        'm_cost': range2rect(size, num_var),
        'power_cost': [x*0 for x in range(size)],                           # 记录每个chrom的能耗代价
        'balance_cost': [x*0 for x in range(size)],                         # 记录每个chrom的负载均衡指数
        'migration_time': [x*0+time_base*num_var for x in range(size)],     # 记录迁移时间，这里指定为固定值
        'rank': [x*0 for x in range(size)],                                 # 记录每个chrom排名，rank值越大，排名越靠后
        'elite_power': 999999.0*num_var,                                    # 记录每代种群中最优秀解的能耗代价值
        'elite_balance': 999999.0*num_var,                                  # 记录每代种群中最优秀解的负载均衡指数
        'elite_migration_time': time_base*num_var,                          # ........的迁移时间
        'elite_chrom': range(num_var)                                       # 保存每代种群中精英chrom
    }
    return population

def init_population(chr1):
    '''initialize chroms in population'''
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            chr1['population'][i][j] = 0
            chr1['p_cost'][i][j] = 0
            chr1['m_cost'][i][j] = 0
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            while(True):
                tmp_position = random.randint(0, chr1['num_var']-1)                # 为第i个chrom中第j个vm随机生成放置pm编号
                tmp_p_cost = chr1['p_cost'][i][tmp_position] + chr1['rp'][j]       # 并统计该标号pm的cpu使用量
                tmp_m_cost = chr1['m_cost'][i][tmp_position] + chr1['rm'][j]       # 统计该标号pm的mem使用量
                if(tmp_p_cost < 1.0 and tmp_m_cost < 1.0):
                    chr1['population'][i][j] = tmp_position
                    chr1['p_cost'][i][tmp_position] = tmp_p_cost
                    chr1['m_cost'][i][tmp_position] = tmp_m_cost
                    break
    for i in range(chr1['size']):
        chr1['population0'][i] = tuple(chr1['population'][i])                      # 对初始化后的种群进行保存
    return chr1

def mbbode_migration(chr1):
    '''migrate SIVs among in population'''
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            rand_sum = random.random()
            lambda_scale = (chr1['lambdaa'][i] - chr1['lambdaa'][0]) / (chr1['lambdaa'][chr1['size']-1] - chr1['lambdaa'][0])   # 标准化迁入率lambda_scale，即该chrom被选为迁入chrom的概率
            if(rand_sum < lambda_scale):                                            # 被选为迁入chrom中的第j元素有rand_sum概率被选为迁出chrom随机生成的SIV替换
                index1 = random.randint(0, chr1['size']-1)
                index2 = random.randint(0, chr1['size']-1)
                while(index1 == index2):
                    index2 = random.randint(0, chr1['size']-1)
                #实现差分迁移，计算被迁入率选中的chrom中，随机选中的SIV将被差分迁移引入的新SIV值（该位置vm将被重新安排的pm标号）
                tmp_position = int(chr1['population'][i][j] + chr1['f']*(chr1['population'][i][j] - chr1['population'][index1][j]) + chr1['f']*(chr1['population'][index1][j] - chr1['population'][index2][j] + 0.5)) % chr1['num_var']
                if(tmp_position < 0):
                    tmp_position = -tmp_position
                while(True):
                    tmp_p_cost = chr1['p_cost'][i][tmp_position] + chr1['rp'][j]
                    tmp_m_cost = chr1['m_cost'][i][tmp_position] + chr1['rm'][j]
                    if(tmp_p_cost < 1.0 and tmp_m_cost < 1.0):                  # 若该vm重新安排至pm(tmp_position)后，cpu,mem资源使用量不超标，则进行替换;否则继续循环查找可容纳该vm的新pm
                        original_position = chr1['population'][i][j]
                        chr1['population'][i][j] = tmp_position
                        chr1['p_cost'][i][tmp_position] = tmp_p_cost
                        chr1['m_cost'][i][tmp_position] = tmp_m_cost
                        chr1['p_cost'][i][original_position] -= chr1['rp'][j]
                        chr1['m_cost'][i][original_position] -= chr1['rm'][j]
                        break
                    else:
                        tmp_position = (tmp_position + 1) % chr1['num_var']
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            if(chr1['p_cost'][i][j] < 0.000001):
                chr1['p_cost'][i][j] = 0
            if(chr1['m_cost'][i][j] < 0.000001):
                chr1['m_cost'][i][j] = 0
    return chr1

def mbbode_mutation(chr1):
    '''mutate SIVs in a population'''
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            rand_num = random.random()
            if(rand_num < chr1['p_mutate']):                                 # chromi的第j元素(vmj) 被突变率选中,随机生成新的SIV值，进行替换
                tmp_position = random.randint(0, chr1['num_var'] - 1)
                while(True):
                    tmp_p_cost = chr1['p_cost'][i][tmp_position] + chr1['rp'][j]
                    tmp_m_cost = chr1['m_cost'][i][tmp_position] + chr1['rm'][j]
                    if(tmp_p_cost < 1.0 and tmp_m_cost < 1.0):
                        original_position = chr1['population'][i][j]
                        chr1['population'][i][j] = tmp_position
                        chr1['p_cost'][i][tmp_position] = tmp_p_cost
                        chr1['m_cost'][i][tmp_position] = tmp_m_cost
                        chr1['p_cost'][i][original_position] -= chr1['rp'][j]
                        chr1['m_cost'][i][original_position] -= chr1['rm'][j]
                        break
                    else:
                        tmp_position = (tmp_position + 1) % chr1['num_var']

    # 将 chr1['p_cost'],chr1['m_cost']清零，并依旧迁移进化后的chr1['chrom']，重新计算（为了确保）
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            chr1['p_cost'][i][j] = 0
            chr1['m_cost'][i][j] = 0
    for i in range(chr1['size']):
        for j in range(chr1['num_var']):
            tmp_position = chr1['population'][i][j]
            chr1['p_cost'][i][tmp_position] += chr1['rp'][j]
            chr1['m_cost'][i][tmp_position] += chr1['rm'][j]
    return chr1

def mbbode_cost(chr1):
    '''compute the costs of three objectives in MBBO/DE algorithm'''
    # compute the power cost of chr1
    for i in range(chr1['size']):
        chr1['power_cost'][i] = 0.0
        for j in range(chr1['num_var']):
            x = chr1['p_cost'][i][j]
            if(x > 0.0):
                chr1['power_cost'][i] += (446.7 + 5.28*x - 0.04747*x*x + 0.000334*x*x*x)             # 计算size个chrom的能耗值  能耗与cpu使用率紧密关系
    # compute the balance cost of chr1
    for i in range(chr1['size']):
        chr1['balance_cost'][i] = 0.0
        load_index = range(chr1['num_var'])
        average_load_index = 0.0
        for j in range(chr1['num_var']):
            load_index[j] = 1.0 / (1.0005 - chr1['p_cost'][i][j]) / (1.0005 - chr1['m_cost'][i][j])  # 各chrom中每个vm的负载指数
            average_load_index += load_index[j]                                                      # 各chrom所有vms的总负载指数
        average_load_index /= 50                                                                     # 各chrom的平均负载指数
        for j in range(chr1['num_var']):
            chr1['balance_cost'][i] = chr1['balance_cost'][i] + (load_index[j] - average_load_index)*(load_index[j] - average_load_index)
        chr1['balance_cost'][i] = math.sqrt(chr1['balance_cost'][i] / 50)                            # 各chrom的负载均衡代价值
    # compute the migration_time of chr1
    for i in range(chr1['size']):
        chr1['migration_time'][i] = 0
        for j in range(chr1['num_var']):
            if(chr1['population'][i][j] != chr1['population0'][i][j]):
                chr1['migration_time'][i] += chr1['time_base']
    return chr1

def mbbode_rank(chr1):
    '''compute each chrom rank of population with non-dominated sorting'''
    # update the value of rank of each chrom
    for i in range(chr1['size']):
        chr1['rank'][i] = 0
    for i in range(chr1['size']):
        for j in range(i+1, chr1['size']):                             # 按照非支配解排序non-dominated sorting
            if(chr1['power_cost'][i] < chr1['power_cost'][j]):
                chr1['rank'][j] += 1
            elif(chr1['power_cost'][i] > chr1['power_cost'][j]):
                chr1['rank'][i] += 1
            if(chr1['balance_cost'][i] < chr1['balance_cost'][j]):
                chr1['rank'][j] += 1
            elif(chr1['balance_cost'][i] > chr1['balance_cost'][j]):
                chr1['rank'][i] += 1
            if(chr1['migration_time'][i] < chr1['migration_time'][j]):
                chr1['rank'][j] += 1
            elif(chr1['migration_time'][i] > chr1['migration_time'][j]):
                chr1['rank'][i] += 1
    rank = 999999
    for i in range(chr1['size']):
        if(rank > chr1['rank'][i]):
            rank = chr1['rank'][i]                                    # 查找最优秀解，rank最小
    for i in range(chr1['size']):
        if(rank == chr1['rank'][i]):                                  # 获取该最优秀解size编号
            if(chr1['elite_power'] > chr1['power_cost'][i] and chr1['elite_balance'] > chr1['balance_cost'][i]):    # 若上一代精英解没有当前代最优解优秀，则替换之
                chr1['elite_power'] = chr1['power_cost'][i]
                chr1['elite_balance'] = chr1['balance_cost'][i]
                chr1['elite_migration_time'] = chr1['migration_time'][i]
                for j in range(chr1['num_var']):
                    chr1['elite_chrom'][j] = chr1['population'][i][j]
            else:
                for j in range(chr1['num_var']):                       # 若当前代最优解没有精英解优秀，则精英解保留不变并用其替换当代最优解（即下标为0的chrom）
                    chr1['population'][0][j] = chr1['elite_chrom'][j]
    return chr1

def mbbode_get_best_chr(chr1, chr2):
    '''get the best chrom from RDD'''
    if(chr1.power_cost < chr2.power_cost and chr1.balance_cost < chr2.balance_cost):
        return chr1
    else:
        return chr2

def get_best_cost(chr1, chr2):
    '''get the best cost from best_chrom'''
    best_cost = [0.0, 0.0, 0.0]
    if(chr1['elite_power'] > chr2['elite_power'] and chr1['elite_balance'] > chr2['elite_balance']):
        best_cost[0] = chr2['elite_power']
        best_cost[1] = chr2['elite_balance']
        best_cost[2] = chr2['elite_migration_time']
    else:
        best_cost[0] = chr1['elite_power']
        best_cost[1] = chr1['elite_balance']
        best_cost[2] = chr1['elite_migration_time']
    return best_cost

if __name__ == '__main__':
    # define some parameters of the MBBO/DE algorithm
    generation = 1000
    size = 10
    num_var = 200
    p_mutate = 0.01
    f = 0.6
    rp_u = 0.25
    rm_u = 0.25
    p = 1.0
    time_base = 65

    # initialize the CPU and Main Memory utilizations of virtual machines'
    rp = range(num_var)                                                  # 记录每个虚拟机对cpu,mem使用情况
    rm = range(num_var)
    for i in range(num_var):
        rp[i] = random.uniform(0, 0.5)
        rm[i] = random.uniform(0, 0.25)
    for i in range(num_var):
        r = random.random()
        if((r < p and rp[i] >= rp_u) or (r >= p and rp[i] < rp_u)):     # 控制一定概率出现的高cpu或高mem情况
            rm[i] += rm_u
        if(rm[i] >= 1.0):
            rm[i] -= 1.0

    # compute the Cosine migration rates
    lambdaa = range(size)
    mu = range(size)
    for i in range(size):
        lambdaa[i] = math.cos(float(size - (i + 1)) / size)             # 迁入率，采用余弦迁移模型，被选中作为迁入岛屿的概率=改岛屿的排名/子population中总chrom个数再求余弦值
        mu[i] = math.sin(float(size - (i + 1)) / size)                  # 迁出率


    # deployed in Spark, here father_size is the number of populations in MBBO/DE
    father_size = 5
    init_population = range(father_size)                                # 初代father_size个populations，each population has num size of chrom候选解
    for i in range(father_size):                                        # 串行初始化5个populations
        init_population[i] = make_population(size, num_var, rp, rm, f, p_mutate, time_base, lambdaa)
        init_population[i] = init_population(init_population[i])
        init_population[i] = mbbode_cost(init_population[i])
        init_population[i] = mbbode_rank(init_population[i])

    sc = SparkContext(appName="Paralleled MBBO/DE algorithm")
    # SparkContext():Main entry point for Spark functionality. A SparkContext represents the connection to a Spark cluster,and can be used to create RDD and broadcast variables on that cluster

    elite_cost = [9999.9*num_var, 9999.9*num_var, time_base*num_var]   # 初设的全局精英解能耗代价、负载均衡指数、迁移时间
    time1 = time.time()
    for g in range(generation):                                        # 设置最大迭代次数
        pop = sc.parallelize(init_population)                          #parallelize()：Distribute a local Python collection to form an RDD
        init_population = pop.map(mbbode_migration).map(mbbode_mutation).map(mbbode_cost).map(mbbode_rank).collect()  #RDD.map():Return a new RDD by applying a function to each element of this RDD
        for i in range(father_size):                                   # 获取全局最优解的能耗代价、负载均衡指数、以及迁移时间
            if(elite_cost[0] > init_population[i]['elite_power'] and elite_cost[1] > init_population[i]['elite_balance']):
                elite_cost[0] = init_population[i]['elite_power']
                elite_cost[1] = init_population[i]['elite_balance']
                elite_cost[2] = init_population[i]['elite_migration_time']
        print elite_cost[0], elite_cost[1], elite_cost[2]

    time2 = time.time()
    print 'Time cost: ', (time2 - time1), '\n'
    sc.stop()
