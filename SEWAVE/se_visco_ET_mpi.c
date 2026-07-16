#include "visco_wave_2d.h"
#include "mpi.h"

int main(int argc, char **argv)
{
    initargs(argc, argv);
    int i, j, k;
    int is, ix, it, iz;

    int flag_homo = 0;              // 1表示取直达波，1表示不取直达波，
    int flag_generate_seis = 0;     // 1表示生成地震记录，0表示不生成地震记录
    int type_compute_Laplace = 1;           // 1表示伪谱法计算拉普拉斯，0表示差分计算拉普拉斯
    int flag_smooth = 1;            // 1表示光滑模型，0表示不光滑模型
    int amp_compensation_sign = -1; // 振幅补偿符号, 1表示衰减，-1表示补偿
    char *PS_or_CS;  // 伪谱法或CS延拓计算分数阶拉普拉斯算子

    char *input_seis_record,*output_seis_record;
    char *mig_result,*mig_temp;

    visco_wave_2d_t p;
    visco_wave_2d_getpar(&p, argc, argv);
    if (!getparstring("PS_or_CS", &PS_or_CS))
        PS_or_CS = "PS";
    if (!getparint("flag_homo", &flag_homo))
        flag_homo = 0;
    if (!getparint("flag_generate_seis", &flag_generate_seis))
        flag_generate_seis = 0;
    if (!getparint("type_compute_Laplace", &type_compute_Laplace))
        type_compute_Laplace = 1;
    if (!getparint("flag_smooth", &flag_smooth))
        flag_smooth = 1;
    if (!getparint("amp_compensation_sign", &amp_compensation_sign))
        amp_compensation_sign = -1;


if (!getparstring("input_seis_record", &input_seis_record))
        input_seis_record = "input.seis";
    if (!getparstring("output_seis_record", &output_seis_record))
        output_seis_record = "output.seis";
    if (!getparstring("mig_result", &mig_result))
        mig_result = "mig.res";
    if (!getparstring("mig_temp", &mig_temp))
        mig_temp = "./temp/mig.temp";

    visco_wave_2d_init(&p, flag_homo, flag_smooth);

    

    // 打印参数
    printf("----------parameters-----------\n");
    printf("nx = %d\n", p.nx);
    printf("nz = %d\n", p.nz);
    printf("dxz = %f\n", p.dxz);
    printf("nt = %d\n", p.nt); // 打印nt
    printf("dt = %f\n", p.dt); // 打印nt
    printf("ns = %d\n", p.ns);
    printf("fs = %f\n", p.fs);
    printf("ds = %f\n", p.ds);
    printf("sx = %f\n", p.sx); // 打印sx
    printf("sz = %f\n", p.sz); // 打印sz
    printf("fdom = %f\n", p.fdom); // 打印f0
    printf("f0 = %f\n", p.f0); // 打印f0
    printf("nbc = %d\n", p.nbc);
    printf("L = %d\n", p.L);
    printf("alpha = %f\n", p.alpha);
    printf("amp_compensation_sign = %d\n",amp_compensation_sign);
    printf("-------------------------------\n");

    int np, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    p.sz = p.sz + p.nbc; // 假设震源深度为nbc
    p.rz = p.sz;         // 假设检波器深度与震源深度相同

    // 注意生成地震记录时不要光滑模型
    if (flag_generate_seis == 1)
    {

        // 写一个炮循环, 生成多炮正演地震记录
        for (is = myid; is < p.ns; is = is + np)
        {
            SE_MESSAGE("Seismic generating Start!!!!! \n");
            SE_MESSAGE("is = : %d \n", is);

            // 观测系统设置
            p.sx = p.fs + is * p.ds + p.nbc;
            printf("p.sx = %f\n", p.sx); // 打印p.sx

            char fname[256];
            sprintf(fname, "./temp/temp_%d", is);
            FILE *fp_allshots = fopen(fname, "wb");

            if (fp_allshots == NULL)
            {
                err("error in opening file");
            }
            // 震源波场正传
            get_ricker_source(&p, p.fdom);

            // 传播
            printf("\n");
            printf("\n");
            printf("Using %s to compute fractional Laplace operator \n", PS_or_CS);
            printf("\n");
            printf("\n");
            
            visco_wave_2d_seismic_generate_TVM_PSandCS(&p, type_compute_Laplace, PS_or_CS);
            


            for (i = 0; i < p.nx - 2 * p.nbc; i++)
            {
                for (j = 0; j < p.nt; j++)
                {
                    fwrite(&p.record[i][j], sizeof(float), 1, fp_allshots);
                }
            }

            fclose(fp_allshots);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (myid == 0)
        {
            // 将temp文件夹内的炮数据按照顺序保存为一个文件
            FILE *fp_allshots = fopen(output_seis_record, "wb");
            if (fp_allshots == NULL)
            {
                err("error in opening file seis_allshots.bin");
            }

            for (is = 0; is < p.ns; is++)
            {
                char fname[256];
                sprintf(fname, "./temp/temp_%d", is);
                FILE *fp_temp = fopen(fname, "rb");
                if (fp_temp == NULL)
                {
                    err("error in opening file %s", fname);
                }

                float *data = (float *)malloc(sizeof(float) * (p.nx - 2 * p.nbc) * p.nt);
                fread(data, sizeof(float), (p.nx - 2 * p.nbc) * p.nt, fp_temp);
                fwrite(data, sizeof(float), (p.nx - 2 * p.nbc) * p.nt, fp_allshots);

                fclose(fp_temp);
                free(data);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        
        // 成像
        float **seis_allshots;
        seis_allshots = alloc2float(p.nt, (p.nx - 2 * p.nbc) * p.ns);
        read_binary_file_2d(input_seis_record, seis_allshots, (p.nx - 2 * p.nbc) * p.ns, p.nt);

        float **receiver_seis;
        receiver_seis = alloc2float(p.nt, p.nx - 2 * p.nbc);
        memset(receiver_seis[0], 0, sizeof(float) * (p.nx - 2 * p.nbc) * p.nt); // 置0初始化

        float **Image, **Fimage, **Fimage_lap;
        Image = alloc2float(p.nz - 2 * p.nbc, p.nx - 2 * p.nbc);
        Fimage = alloc2float(p.nz - 2 * p.nbc, p.nx - 2 * p.nbc);
        memset(Image[0], 0, sizeof(float) * (p.nx - 2 * p.nbc) * (p.nz - 2 * p.nbc));
        memset(Fimage[0], 0, sizeof(float) * (p.nx - 2 * p.nbc) * (p.nz - 2 * p.nbc));
        Fimage_lap = alloc2float(p.nz - 2 * p.nbc, p.nx - 2 * p.nbc);
        memset(Fimage_lap[0], 0, sizeof(float) * (p.nx - 2 * p.nbc) * (p.nz - 2 * p.nbc));
        float **data;
        for (is = myid; is < p.ns; is = is + np)
        {
            SE_MESSAGE("Imageing Start!!!!! \n");
            SE_MESSAGE("is = : %d \n", is);

            // 初始化波场快照
            float ***forward_wavefield;
            forward_wavefield = alloc3float(p.nt, p.nz, p.nx);
            memset(forward_wavefield[0][0], 0, sizeof(float) * p.nt * p.nx * p.nz);

            float **Image_1shot;
            Image_1shot = alloc2float(p.nz - 2 * p.nbc, p.nx - 2 * p.nbc);
            memset(Image_1shot[0], 0, sizeof(float) * (p.nx - 2 * p.nbc) * (p.nz - 2 * p.nbc));

            // 观测系统设置
            p.sx = p.fs + is * p.ds + p.nbc;
            printf("p.sx = %f\n", p.sx); // 打印p.sx

            // 传播
            printf("\n");
            printf("\n");
            printf("Using %s to compute fractional Laplace operator \n", PS_or_CS);
            printf("\n");
            printf("\n");
            // 震源波场正传
            get_ricker_source(&p, p.fdom);
            visco_wave_2d_propagation_TVM_PSandCS(&p, 0, type_compute_Laplace, PS_or_CS, amp_compensation_sign, forward_wavefield, Image_1shot); 

            // 获取每炮地震记录
            for (i = 0; i < p.nx - 2 * p.nbc; i++)
            {
                for (j = 0; j < p.nt; j++)
                {
                    receiver_seis[i][j] = seis_allshots[i + is * (p.nx - 2 * p.nbc)][j];
                }
            }
            // write_binary_file_2d("receiver_seis111.bin", receiver_seis, p.nx - 2 * p.nbc, p.nt);
            // exit(0);

            // 检波器波场反传
            get_receiver_seis(&p, receiver_seis);
            visco_wave_2d_propagation_TVM_PSandCS(&p, 1, type_compute_Laplace, PS_or_CS, amp_compensation_sign, forward_wavefield, Image_1shot);

            // 多炮叠加
            for (ix = 0; ix < p.nx - 2 * p.nbc; ix++)
            {
                for (iz = 0; iz < p.nz - 2 * p.nbc; iz++)
                {
                    Image[ix][iz] = Image[ix][iz] + Image_1shot[ix][iz];
                }
            }


            char fname[256];
            sprintf(fname, "%s_%d",mig_temp, is);
            FILE *fp_temp1 = fopen(fname, "wb");
            if (fp_temp1 == NULL)
            {
                err("error in opening file %s", fname);
            }

            for (ix = 0; ix < p.nx - 2 * p.nbc; ix++)
            {
                for (iz = 0; iz < p.nz - 2 * p.nbc; iz++)
                {
                    fwrite(&Image_1shot[ix][iz], sizeof(float), 1, fp_temp1);
                }
            }

            free3float(forward_wavefield);
            free2float(Image_1shot);
            fclose(fp_temp1);
        }

        MPI_Reduce(Image[0], Fimage[0], (p.nx - 2 * p.nbc) * (p.nz - 2 * p.nbc), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);


        if (myid == 0)
        {
            
            // 输出互相关成像结果
            FILE *fp_image;
            fp_image = fopen(mig_result, "wb");
            if (fp_image == NULL)
            {
                err("error in opening Image.bin");
            }
            for (ix = 0; ix < p.nx - 2 * p.nbc; ix++)
            {
                for (iz = 0; iz < p.nz - 2 * p.nbc; iz++)
                {
                    fwrite(&Fimage[ix][iz], sizeof(float), 1, fp_image);
                }
            }
            fclose(fp_image);

                        
            // // 输出拉普拉斯滤波成像结果
            // read_binary_file_2d("./result/mar/Image_mar_10shot_PS.bin", Image, p.nx - 2 * p.nbc, p.nz - 2 * p.nbc);
            // laplace_denoise(Image, Fimage_lap, p.nx - 2 * p.nbc, p.nz - 2 * p.nbc);
            // write_binary_file_2d("./result/mar/Image_mar_10shot_PS_lap.bin", Fimage_lap, p.nx - 2 * p.nbc, p.nz - 2 * p.nbc);

            SE_MESSAGE("Image saved successful!!\n");


        }

        // 释放内存
        free2float(seis_allshots);
        free2float(receiver_seis);
        free2float(Image);
        free2float(Fimage);
        free2float(Fimage_lap);
        
    }

    visco_wave_2d_free(&p);

    MPI_Finalize();

    return 0;
}
