#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
Created on 2016.9.20
@author: Amy
'''
import time
import random
import math

if __name__ == '__main__':
    time1 = time.time()

    # define some parameters of the MBBO/DE algorithm
    generation = 2000                                                        # num of iteration
    size = 50                                                                # 50 chroms in every population
    num_var = 200                                                            # num of vms need to be arranged
    p_mutate = 0.01                                                          # the probability of the mutate
    f = 0.6                                                                  # difference factor
    rp_u = 0.25                                                              # the used rp
    rm_u = 0.25                                                              # the used rm
    p = 1.0
    time_base = 65

    # initialize the CPU and Main Memory utilizations of virtual machines'
    rp = range(num_var)
    rm = range(num_var)
    for i in range(num_var):
        rp[i] = random.uniform(0, 0.5)
        rm[i] = random.uniform(0, 0.25)
    for i in range(num_var):
        r = random.random()
        if((r < p and rp[i] >= rp_u) or (r >= p and rp[i] < rp_u)):
            rm[i] += rm_u
        if(rm[i] >= 1.0):                                                   # 对随机生成的各vm资源请求进行合理化
            rm[i] -= 1.0

    # initialize the init_population
    init_population = range(size)
    p_cost = range(size)
    m_cost = range(size)
    for i in range(size):
        init_population[i] = range(num_var)                                 # init_population包含50个list，每个list有200个元素
        p_cost[i] = range(num_var)
        m_cost[i] = range(num_var)
    for i in range(size):
        for j in range(num_var):
            init_population[i][j] = 0
            p_cost[i][j] = 0
            m_cost[i][j] = 0
    for i in range(size):
        for j in range(num_var):
            while(True):
                tmp_position = random.randint(0, num_var-1)                  # the nums of pm which can be used is the same as nums of vms
                tmp_p_cost = p_cost[i][tmp_position] + rp[j]
                tmp_m_cost = m_cost[i][tmp_position] + rm[j]
                if(tmp_p_cost < 1.0 and tmp_m_cost < 1.0):
                    init_population[i][j] = tmp_position
                    p_cost[i][tmp_position] = tmp_p_cost                     # p_cost is the used-cpu each chrom in population
                    m_cost[i][tmp_position] = tmp_m_cost                     # p_mem is the mem
                    break

    # define another population to be used in the process of optimization,
    # the init_population is used as a comprasion.
    tmp_population = range(size)
    for i in range(size):
        tmp_population[i] = range(num_var)
    for i in range(size):
        for j in range(num_var):
            tmp_population[i][j] = init_population[i][j]                    # save inital 50x200chrom

    # compute the power cost
    power_cost = range(size)
    for i in range(size):
        power_cost[i] = 0.0
    for i in range(size):
        for j in range(num_var):
            tmp_p_cost = p_cost[i][j]
            if(tmp_p_cost > 0.0):
                power_cost[i] += (446.7 + 5.28*tmp_p_cost - 0.04747*tmp_p_cost*tmp_p_cost + 0.000334*tmp_p_cost*tmp_p_cost*tmp_p_cost)  #compute power_cost of each chrom in population

    # compute the balance cost
    balance_cost = range(size)
    load_index = range(size)
    count_of_active_machine = range(size)
    average_load_index = range(size)
    for i in range(size):
        balance_cost[i] = 0
        count_of_active_machine[i] = 0
        average_load_index[i] = 0
        load_index[i] = range(num_var)
    for i in range(size):
        for j in range(num_var):
            load_index[i][j] = 1.0 / (1.0005 - p_cost[i][j]) / (1.0005 - m_cost[i][j])
            average_load_index[i] += load_index[i][j]
        average_load_index[i] /= size
    for i in range(size):
        for j in range(num_var):
            balance_cost[i] = balance_cost[i] + (load_index[i][j] - average_load_index[i])*(load_index[i][j] - average_load_index[i])
        balance_cost[i] = math.sqrt(balance_cost[i] / size)

    # compute the migration time
    time_of_migration = range(size)
    for i in range(size):
        time_of_migration[i] = time_base*num_var   #我靠，迁移时间设置的固定值啊


    # sort the population by NDRS（non-dominated sorting） method
    #非支配解排序：设任何二解S1 及S2 对所有目标而言，S1均优于S2，则我们称S1 支配S2，若S1 的解没有被其他解所支配，则S1 称为非支配解
    #为解决基于帕累托(Pareto)支配解排序的多目标进化算法高时间复杂度问题,依据非支配解排序潜在特性,介绍了一种快速的非支配解排序方法,
    #每次只处理当前种群中最高等级个体,且在分配等级的同时,能选择个体进入下一代,下一代被选足时即结束程序,减少了排序处理个体的数量,大幅度降低时间复杂度
    rank_of_population = range(size)
    for i in range(size):
        rank_of_population[i] = 0
    for i in range(size):
        for j in range(size):
            if(power_cost[i] < power_cost[j]):
                rank_of_population[j] += 1
            elif(power_cost[i] > power_cost[j]):
                rank_of_population[i] += 1
            if(balance_cost[i] < balance_cost[j]):
                rank_of_population[j] += 1
            elif(balance_cost[i] > balance_cost[j]):           #代价值越大，rank值越高，排名越靠后
                rank_of_population[i] += 1
            if(time_of_migration[i] < time_of_migration[j]):
                rank_of_population[j] += 1
            elif(time_of_migration[i] > time_of_migration[j]):
                rank_of_population[i] += 1

    # find & save the init-solution, init-cost, and the best-solution
    best_population = range(num_var)
    tmp_best_population = range(num_var)
    min_rank = 500
    for i in range(size):
        if(min_rank > rank_of_population[i]):
            min_rank = rank_of_population[i]
    for i in range(size):
        if(min_rank == rank_of_population[i]):
            init_best_power_cost = power_cost[i]
            init_best_balance_cost = balance_cost[i]
            init_best_migration_time = time_of_migration[i]
            best_power_cost = power_cost[i]
            best_balance_cost = balance_cost[i]
            best_migration_time = time_of_migration[i]
            for j in range(num_var):
                best_population[j] = tmp_population[i][j]                  # save the best chrom in inital population which has the smallest rank meaning the lowest cost，and save in best_population
            break

    #-------------------------------------------------------------------
    #--------------------start--optimization----------------------------开始进行迁移、突变等进化
    # compute migration rates
    lambdaa = range(size)
    mu = range(size)
    for i in range(size):
        lambdaa[i] = math.cos(float(size - (i + 1)) / size)               # use Cosine model to compute immigration probability：lambdaa，and emmigration probability:mu
        mu[i] = math.sin(float(size - (i + 1)) / size)

    for g in range(generation):
        # migration SIVs among populations with difference migration
        for i in range(size):
            for j in range(num_var):
                rand_sum = random.random()
                lambda_scale = (lambdaa[i] - lambdaa[0]) / (lambdaa[size - 1] - lambdaa[0])   #标准化迁入率lambdaa
                if(rand_sum < lambda_scale):                             # each chrom has a probability to be choosed replace its j-th SIV
                    index1 = random.randint(0, size-1)
                    index2 = random.randint(0, size-1)
                    while(index1 == index2):
                        index2 = random.randint(0, size-1)
                    tmp_position = int(tmp_population[i][j] + f*(tmp_population[i][j] - tmp_population[index1][j]) + f*(tmp_population[index1][j] - tmp_population[index2][j] + 0.5)) % num_var    #实现差分迁移，计算被选为迁入chrom的迁入SIV的位置
                    if(tmp_position < 0):
                        tmp_position = -tmp_position
                    while(True):
                        tmp_p_cost = p_cost[i][tmp_position] + rp[j]
                        tmp_m_cost = m_cost[i][tmp_position] + rm[j]
                        if(tmp_p_cost < 1.0 and tmp_m_cost < 1.0):
                            original_position = tmp_population[i][j]       # tmp_population is initail 50x200 population
                            tmp_population[i][j] = tmp_position
                            p_cost[i][tmp_position] = tmp_p_cost
                            m_cost[i][tmp_position] = tmp_m_cost
                            p_cost[i][original_position] -= rp[j]
                            m_cost[i][original_position] -= rm[j]
                            break
                        else:
                            tmp_position = (tmp_position + 1) % num_var

        # recompute the boundaries of p_cost and m_cost
        for i in range(size):
            for j in range(num_var):
                if(p_cost[i][j] < 0):
                    p_cost[i][j] = 0
                if(m_cost[i][j] < 0):
                    m_cost[i][j] = 0

        # mutate SIVs among populations
        for i in range(size):
            for j in range(num_var):
                rand_num = random.random()
                if(rand_num < p_mutate):
                    tmp_position = random.randint(0, num_var-1)
                    while(True):
                        tmp_p_cost = p_cost[i][tmp_position] + rp[j]
                        tmp_m_cost = m_cost[i][tmp_position] + rm[j]
                        if(tmp_p_cost < 1.0 and tmp_m_cost < 1.0):
                            original_position = tmp_population[i][j]
                            tmp_population[i][j] = tmp_position
                            p_cost[i][tmp_position] = tmp_p_cost
                            m_cost[i][tmp_position] = tmp_m_cost
                            p_cost[i][original_position] -= rp[j]
                            m_cost[i][original_position] -= rm[j]
                            break
                        else:
                            tmp_position = (tmp_position + 1) % num_var

        # recompute the p_cost and m_cost
        for i in range(size):
            for j in range(num_var):
                p_cost[i][j] = 0
                m_cost[i][j] = 0
        for i in range(size):
            for j in range(num_var):
                tmp_position = tmp_population[i][j]
                p_cost[i][tmp_position] += rp[j]
                m_cost[i][tmp_position] += rm[j]

        # recompute the power cost
        for i in range(size):
            power_cost[i] = 0
        for i in range(size):
            for j in range(num_var):
                tmp_p_cost = p_cost[i][j]
                if(tmp_p_cost > 0):
                    power_cost[i] += (446.7 + 5.28*tmp_p_cost - 0.04747*tmp_p_cost*tmp_p_cost + 0.000334*tmp_p_cost*tmp_p_cost*tmp_p_cost)

        # recompute the balance cost
        for i in range(size):
            balance_cost[i] = 0
            count_of_active_machine[i] = 0
            average_load_index[i] = 0
            for j in range(num_var):
                load_index[i][j] = 0
        for i in range(size):
            for j in range(num_var):
                load_index[i][j] = 1.0 / (1.0005 - p_cost[i][j]) / (1.0005 - m_cost[i][j])
                average_load_index[i] += load_index[i][j]
            average_load_index[i] /= size
        for i in range(size):
            for j in range(num_var):
                balance_cost[i] = balance_cost[i] + (load_index[i][j] - average_load_index[i])*(load_index[i][j] - average_load_index[i])
            balance_cost[i] = math.sqrt(balance_cost[i] / size)

        # recompute the migration time
        for i in range(size):
            time_of_migration[i] = 0
            for j in range(num_var):
                if(tmp_population[i][j] != init_population[i][j]):
                    time_of_migration[i] += 65

        # sort the new populations
        for i in range(size):
            rank_of_population[i] = 0
        for i in range(size):
            for j in range(size):
                if(power_cost[i] < power_cost[j]):
                    rank_of_population[j] += 1
                elif(power_cost[i] > power_cost[j]):
                    rank_of_population[i] += 1
                if(balance_cost[i] < balance_cost[j]):
                    rank_of_population[j] += 1
                elif(balance_cost[i] > balance_cost[j]):
                    rank_of_population[i] += 1
                if(time_of_migration[i] < time_of_migration[j]):
                    rank_of_population[j] += 1
                elif(time_of_migration[i] > time_of_migration[j]):
                    rank_of_population[i] += 1

        # find the new best population
        min_rank = 500
        for i in range(size):
            if(min_rank > rank_of_population[i]):
                min_rank = rank_of_population[i]
        for i in range(size):
            if(min_rank == rank_of_population[i]):
                tmp_best_power_cost = power_cost[i]
                tmp_best_balance_cost = balance_cost[i]
                tmp_best_migration_time = time_of_migration[i]
                tmp_best_index = i
                break

        # determine whether update the best population or not
        if(tmp_best_power_cost < best_power_cost and tmp_best_balance_cost < best_balance_cost):
            best_power_cost = tmp_best_power_cost
            best_balance_cost = tmp_best_balance_cost
            best_migration_time = tmp_best_migration_time
            for j in range(num_var):
                best_population[j] = tmp_population[tmp_best_index][j]
        else:
            for j in range(num_var):
                tmp_population[0][j] = best_population[j]

        # print the process of optimization
        print 'Generation:', g
        print 'Init - Power: ', init_best_power_cost, '; Balance: ', init_best_balance_cost, '; Time: ', init_best_migration_time
        print 'Last - Power: ', best_power_cost, '; Balance: ', best_balance_cost, '; Time: ', best_migration_time
        print 'Curr - Power: ', tmp_best_power_cost, '; Balance: ', tmp_best_balance_cost, '; Time: ', tmp_best_migration_time

    time2 = time.time()
    print 'Time cost: ', (time2 - time1), '\n'
