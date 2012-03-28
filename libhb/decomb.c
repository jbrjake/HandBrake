/* $Id: decomb.c,v 1.14 2008/04/25 5:00:00 jbrjake Exp $

   This file is part of the HandBrake source code.
   Homepage: <http://handbrake.fr/>.
   It may be used under the terms of the GNU General Public License. 
   
   The yadif algorithm was created by Michael Niedermayer.
   Tritical's work inspired much of the comb detection code:
   http://web.missouri.edu/~kes25c/
*/

/**
@file
 - Parameters:
    - Mode : Spatial metric : Motion thresh : Spatial thresh : Mask Filter Mode :
    - Block thresh : Block width : Block height

 - Appended for EEDI2:
    - Magnitude thresh : Variance thresh : Laplacian thresh : Dilation thresh :
    - Erosion thresh : Noise thresh : Max search distance : Post-processing

 - Plus:
    - Parity
    
 - Defaults:
    - 7:2:6:9:1:80:16:16:10:20:20:4:2:50:24:1:-1

Modes can be layered. For example, Yadif (1) + EEDI2 (8) = 9,
which will feed EEDI2 interpolations to yadif.

- Working combos: 
  -  1: Just yadif
  -  2: Just blend
  -  3: Switch between yadif and blend
  -  4: Just cubic interpolate
  -  5: Cubic->yadif
  -  6: Switch between cubic and blend
  -  7: Switch between cubic->yadif and blend
  -  8: Just EEDI2 interpolate
  -  9: EEDI2->yadif
  - 10: Switch between EEDI2 and blend
  - 11: Switch between EEDI2->yadif and blend
  - 17: Yadif->mcdeint
  - 18: Blend->mcdeint
  - 19: Switch between blending and yadif -> mcdeint
  - 20: Cubic->mdeint
  - 21: Cubic->yadif->mcdeint
  - 22: Cubic or blend -> mcdeint
  - 23: Cubic->yadif or blend -> mcdeint
  - 24: EEDI2->mcdeint
  - 25: EEDI2->yadif->mcdeint
  - ...okay I'm getting bored now listing all these different modes
  - 32: Passes through the combing mask for every combed frame (white for combed pixels, otherwise black)
  - 33+: Overlay the combing mask for every combed frame on top of the filtered output (white for combed pixels)
  - 12-15: EEDI2 will override cubic interpolation
  - 16: DOES NOT WORK BY ITSELF-- mcdeint needs to be fed by another deinterlacer
*/

#define MODE_YADIF       1 ///< Use yadif
#define MODE_BLEND       2 ///< Use blending interpolation
#define MODE_CUBIC       4 ///< Use cubic interpolation
#define MODE_EEDI2       8 ///< Use EEDI2 interpolation
#define MODE_MCDEINT    16 ///< Post-process with mcdeint
#define MODE_MASK       32 ///< Output combing masks instead of pictures
#define MODE_GAMMA      128 ///< Scale gamma when decombing
#define MODE_FILTER     256 ///< Filter combing mask
#define MODE_COMPOSITE  512 ///< Overlay combing mask onto picture

#define FILTER_CLASSIC 1 ///< Uses the old way of filtering noise from the combing mask
#define FILTER_ERODE_DILATE 2 ///< Calls the new erode->dilate->filter workflow to remove noise from the mask

#include "hb.h"
#include "hbffmpeg.h"
#include "mpeg2dec/mpeg2.h"
#include "eedi2.h"

#define SUPPRESS_AV_LOG

#define PARITY_DEFAULT   -1

#define MCDEINT_MODE_DEFAULT   -1
#define MCDEINT_QP_DEFAULT      1

#define ABS(a) ((a) > 0 ? (a) : (-(a)))
#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

// Some names to correspond to the pv->eedi_half array's contents
#define SRCPF 0
#define MSKPF 1
#define TMPPF 2
#define DSTPF 3
// Some names to correspond to the pv->eedi_full array's contents
#define DST2PF 0
#define TMP2PF2 1
#define MSK2PF 2
#define TMP2PF 3
#define DST2MPF 4

struct yadif_arguments_s {
    uint8_t **dst;
    int parity;
    int tff;
    int stop;
    int is_combed;
};

struct decomb_arguments_s {
    int stop;
};

struct eedi2_arguments_s {
    int stop;
};

typedef struct yadif_arguments_s yadif_arguments_t;
typedef struct decomb_arguments_s decomb_arguments_t;
typedef struct eedi2_arguments_s eedi2_arguments_t;

typedef struct eedi2_thread_arg_s {
    hb_filter_private_t *pv;
    int plane;
} eedi2_thread_arg_t;

typedef struct decomb_thread_arg_s {
    hb_filter_private_t *pv;
    int segment;
} decomb_thread_arg_t;

typedef struct yadif_thread_arg_s {
    hb_filter_private_t *pv;
    int segment;
} yadif_thread_arg_t;

struct hb_filter_private_s
{
    int              pix_fmt;
    int              width[3];
    int              height[3];

    // Decomb parameters
    int              mode;
    int              filter_mode;
    int              spatial_metric;
    int              motion_threshold;
    int              spatial_threshold;
    int              block_threshold;
    int              block_width;
    int              block_height;
    
    float            gamma_lut[256];
    
    // EEDI2 parameters
    int              magnitude_threshold;
    int              variance_threshold;
    int              laplacian_threshold;
    int              dilation_threshold;
    int              erosion_threshold;
    int              noise_threshold;
    int              maximum_search_distance;
    int              post_processing;

    int              parity;
    int              tff;
    
    int              yadif_ready;

    int              mcdeint_mode;
    int              mcdeint_qp;

    int              mcdeint_outbuf_size;
    uint8_t        * mcdeint_outbuf;
    AVCodecContext * mcdeint_avctx_enc;
    AVFrame        * mcdeint_frame;
    AVFrame        * mcdeint_frame_dec;

    int              deinterlaced_frames;
    int              blended_frames;
    int              unfiltered_frames;

    uint8_t        * ref[4][3];
    int              ref_stride[3];

    /* Make buffers to store comb masks. */
    uint8_t        * mask[3];
    uint8_t        * mask_filtered[3];
    uint8_t        * mask_temp[3];
    int              mask_box_x;
    int              mask_box_y;
    uint8_t          mask_box_color;


    uint8_t        * eedi_half[4][3];
    uint8_t        * eedi_full[5][3];
    int            * cx2;
    int            * cy2;
    int            * cxy;
    int            * tmpc;
    
    AVPicture        pic_in;
    AVPicture        pic_out;
    hb_buffer_t *    buf_out[2];
    hb_buffer_t *    buf_settings;
    
    int              cpu_count;

    hb_thread_t    ** yadif_threads;         // Threads for Yadif - one per CPU
    hb_lock_t      ** yadif_begin_lock;      // Thread has work
    hb_lock_t      ** yadif_complete_lock;   // Thread has completed work
    yadif_arguments_t *yadif_arguments;      // Arguments to thread for work
    
    hb_thread_t    ** decomb_threads;        // Threads for comb detection - one per CPU
    hb_lock_t      ** decomb_begin_lock;     // Thread has work
    hb_lock_t      ** decomb_complete_lock;  // Thread has completed work
    decomb_arguments_t *decomb_arguments;    // Arguments to thread for work

    hb_thread_t    ** eedi2_threads;        // Threads for eedi2 - one per plane
    hb_lock_t      ** eedi2_begin_lock;     // Thread has work
    hb_lock_t      ** eedi2_complete_lock;  // Thread has completed work
    eedi2_arguments_t *eedi2_arguments;    // Arguments to thread for work

//    int              alternator;           // for bobbing parity when framedoubling
};

hb_filter_private_t * hb_decomb_init( int pix_fmt,
                                           int width,
                                           int height,
                                           char * settings );

int hb_decomb_work(      const hb_buffer_t * buf_in,
                         hb_buffer_t ** buf_out,
                         int pix_fmt,
                         int width,
                         int height,
                         hb_filter_private_t * pv );

void hb_decomb_close( hb_filter_private_t * pv );

hb_filter_object_t hb_filter_decomb =
{
    FILTER_DECOMB,
    "Decomb",
    NULL,
    hb_decomb_init,
    hb_decomb_work,
    hb_decomb_close,
};

/**
 @brief Guesses the intensity of a pixel based on 4 neighbors.
 
 @param y0 The pixel 2 above the one being interpolated, weighted -1 or -3
 @param y1 The pixel above the one being interpolated, weighted 5 or 23
 @param y2 The pixel below the one being interpolated, weighted 5 or 23
 @param y3 The pixel 2 below the one being interpolated, weighted -1 or -3
 
 @returns An interpolated value for a pixel, clamped to uint8_t
 
 The weights vary because the first values are for a fast interpolated, which is divided cleanly by 8, and the second values are for the precise interpolation, which is divided by a messy 40.
*/
int cubic_interpolate_pixel( int y0, int y1, int y2, int y3 )
{
    /** Fast but less precise
    int result = ( y0 * -1 ) + ( y1 * 5 ) + ( y2 * 5 ) + ( y3 * -1 );
    result /= 8;
    */
    
    int result = ( y0 * -3 ) + ( y1 * 23 ) + ( y2 * 23 ) + ( y3 * -3 );
    result /= 40;
    
    if( result > 255 )
    {
        result = 255;
    }
    else if( result < 0 )
    {
        result = 0;
    }
    
    return result;
}

/**
 @brief Takes a line of pixels, 
*/
static void cubic_interpolate_line( uint8_t *dst,
                               uint8_t *cur,
                               int plane,
                               int y,
                               hb_filter_private_t * pv )
{
    int w = pv->width[plane];
    int refs = pv->ref_stride[plane];
    int x;

    for( x = 0; x < w; x++)
    {
        int a, b, c, d;
        a = b = c = d = 0;
        
        if( y >= 3 )
        {
            /* Normal top*/
            a = cur[-3*refs];
            b = cur[-refs];
        }
        else if( y == 2 || y == 1 )
        {
            /* There's only one sample above this pixel, use it twice. */
            a = cur[-refs];
            b = cur[-refs];
        }
        else if( y == 0 )
        {
            /* No samples above, triple up on the one below. */
            a = cur[+refs];
            b = cur[+refs];
        }
        
        if( y <= ( pv->height[plane] - 4 ) )
        {
            /* Normal bottom*/
            c = cur[+refs];
            d = cur[3*refs];            
        }
        else if( y == ( pv->height[plane] - 3 ) || y == ( pv->height[plane] - 2 ) )
        {
            /* There's only one sample below, use it twice. */
            c = cur[+refs];
            d = cur[+refs];
        }
        else if( y == pv->height[plane] - 1)
        {
            /* No samples below, triple up on the one above. */
            c = cur[-refs];
            d = cur[-refs];
        }
        
        dst[0] = cubic_interpolate_pixel( a, b, c, d );
        
        dst++;
        cur++;
    }
}

void draw_mask_box( hb_filter_private_t * pv )
{
    int x = pv->mask_box_x;
    int y = pv->mask_box_y;
    int box_width = pv->block_width;
    int box_height = pv->block_height;
    int stride = pv->ref_stride[0];
    uint8_t * mskp = ( pv->mode & MODE_FILTER ) ? pv->mask_filtered[0] : pv->mask[0];

    
    int block_x, block_y;
    for( block_x = 0; block_x < box_width; block_x++)
    {
        mskp[y*stride+x+block_x] = 128;
        mskp[(y+box_height)*stride+x+block_x] = 128;
    }
    
    for( block_y = 0; block_y < box_height; block_y++)
    {
        mskp[stride*(y+block_y)+x] = 128;
        mskp[stride*(y+block_y) + x + box_width] = 128;
    }
    
}

void apply_mask_line( uint8_t * srcp,
                      uint8_t * mskp,
                      int width )
{
    int x;
    
    for( x = 0; x < width; x++ )
    {
        if( mskp[x] == 255 )
        {
            srcp[x] = 255;
        }
        if( mskp[x] == 128 )
        {
            srcp[x] = 128;
        }
    }
}

void apply_mask( hb_filter_private_t * pv )
{
    
    /* draw_boxes */
    draw_mask_box( pv );
    
    int plane, height;
    
    for( plane = 0; plane < 3; plane++ )
    {
        uint8_t * srcp = ( pv->mode & MODE_MCDEINT ) ? pv->pic_in.data[plane] : pv->pic_out.data[plane];
        uint8_t * mskp = ( pv->mode & MODE_FILTER ) ? pv->mask_filtered[plane] : pv->mask[plane];
        
        for( height = 0; height < pv->height[plane]; height++ )
        {
            
            if( !(pv->mode & MODE_COMPOSITE) && plane == 0 )
            {
                memcpy( srcp, mskp, pv->width[plane] );
            }
            else if( !(pv->mode & MODE_COMPOSITE) )
            {
                memset( srcp, 128, pv->width[plane] );
            }
            else if( plane == 0 )
            {
                apply_mask_line( srcp, mskp, pv->width[plane] );
            }
            
//            if( ( pv->mode == MODE_MASK ||
//                    ( pv->mode == ( MODE_MASK + MODE_GAMMA + MODE_FILTER) ) || 
//                    ( pv->mode == ( MODE_MASK + MODE_FILTER ) ) || 
//                    ( pv->mode == ( MODE_MASK + MODE_GAMMA ) ) )
//                  && plane == 0 )
//            {
//                memcpy( srcp, mskp, pv->width[plane] );
//            }
//            else if( pv->mode == MODE_MASK || 
//                   ( pv->mode == ( MODE_MASK + MODE_GAMMA + MODE_FILTER ) ) || 
//                   ( pv->mode == ( MODE_MASK + MODE_FILTER ) ) ||
//                   ( pv->mode == ( MODE_MASK + MODE_GAMMA ) ) )
//            {
//                memset( srcp, 128, pv->width[plane] );
//            }
//            else if( plane == 0 )
//            {
//                apply_mask_line( srcp, mskp, pv->width[plane] );
//            }

            srcp += pv->pic_out.linesize[plane];
            mskp += pv->ref_stride[plane];
        }
    }
    
}

static void store_ref( const uint8_t ** pic,
                             hb_filter_private_t * pv )
{
    memcpy( pv->ref[3],
            pv->ref[0],
            sizeof(uint8_t *)*3 );

    memmove( pv->ref[0],
             pv->ref[1],
             sizeof(uint8_t *)*3*3 );

    int i;
    for( i = 0; i < 3; i++ )
    {
        const uint8_t * src = pic[i];
        uint8_t * ref = pv->ref[2][i];

        int w = pv->width[i];
        int h = pv->height[i];
        int ref_stride = pv->ref_stride[i];

        int y;
        for( y = 0; y < h; y++ )
        {
            memcpy(ref, src, w);
            src = (uint8_t*)src + w;
            ref = (uint8_t*)ref + ref_stride;
        }
    }
}

/* This function may be useful in the future, if we want to output
   a reference to an AVPicture, since they have different strides.
static void get_ref( uint8_t ** pic, hb_filter_private_t * pv, int frm )
{
    int i;
    for( i = 0; i < 3; i++ )
    {
        uint8_t * dst = pic[i];
        const uint8_t * ref = pv->ref[frm][i];
        int w = pv->width[i];
        int ref_stride = pv->ref_stride[i];
        
        int y;
        for( y = 0; y < pv->height[i]; y++ )
        {
            memcpy(dst, ref, w);
            dst += w;
            ref += ref_stride;
        }
    }
}
*/

int blend_filter_pixel( int up2, int up1, int current, int down1, int down2 )
{
    /* Low-pass 5-tap filter */
    int result = 0;
    result += -up2;
    result += up1 * 2;
    result += current * 6;
    result += down1 *2;
    result += -down2;
    result /= 8;

    if( result > 255 )
    {
        result = 255;
    }
    if( result < 0 )
    {
        result = 0;
    }
    
    return result;
}

static void blend_filter_line( uint8_t *dst,
                               uint8_t *cur,
                               int plane,
                               int y,
                               hb_filter_private_t * pv )
{
    int w = pv->width[plane];
    int refs = pv->ref_stride[plane];
    int x;

    for( x = 0; x < w; x++)
    {
        int a, b, c, d, e;
        
        a = cur[-2*refs];
        b = cur[-refs];
        c = cur[0];
        d = cur[+refs];
        e = cur[2*refs];
        
        if( y == 0 )
        {
            /* First line, so A and B don't exist.*/
            a = cur[0];
            b = cur[0];
        }
        else if( y == 1 )
        {
            /* Second line, no A. */
            a = cur[-refs];
        }
        else if( y == (pv->height[plane] - 2) )
        {
            /* Second to last line, no E. */
            e = cur[+refs];
        }
        else if( y == (pv->height[plane] -1) )
        {
            /* Last line, no D or E. */
            d = cur[0];
            e = cur[0];
        }
                
        dst[0] = blend_filter_pixel( a, b, c, d, e );

        dst++;
        cur++;
    }
}

/*
void erode_combing_mask( hb_filter_private_t * pv )
{
    //Take in mask_temp, output to mask_filtered
    int k;
    
    int count;
    
    int erosion_threshold = 2;
    
    int ref;
    for( k = 0; k < 1; k++ )
    {
        ref = pv->ref_stride[k];
        
        int height_limit = pv->height[k] -1;
        int width_limit = pv->width[k] -1;
        
        int y;
        for( y = 1; y < height_limit; ++y )
        {
            uint8_t * cur = &pv->mask_temp[k][y*ref];
            uint8_t * dst = &pv->mask_filtered[k][y*ref];
            
            int x;
            for( x = 1; x < width_limit; ++x )
            {
                dst[0] = 0;
                if( cur[0] == 255 )
                {
                    count = 0;
                    if( cur[-1-ref] == 255 )
                        ++count;
                    if( cur[-ref] == 255 )
                        ++count;
                    if( cur[1-ref] == 255 )
                        ++count;
                    if( cur[-1] == 255 )
                        ++count;
                    if( cur[1] == 255 )
                        ++count;
                    if( cur[-1+ref] == 255 )
                        ++count;
                    if( cur[+ref] == 255 )
                        ++count;
                    if( cur[1+ref] == 255 )
                        ++count;

                    if( count >= erosion_threshold )
                    {
                        dst[0] = 255;
                    }
                }
                cur++;
                dst++;
            }
        }
    }
}

void dilate_combing_mask( hb_filter_private_t * pv )
{
    //Take in mask_temp, output to mask
    int k;
    
    int count;
    
    int dilation_threshold = 4;
    
    int ref;
    for( k = 0; k < 1; k++ )
    {
        ref = pv->ref_stride[k];
        int height_limit = pv->height[k] -1;
        int width_limit = pv->width[k] -1;
        
        int y;
        for( y = 1; y < height_limit; y++ )
        {
            uint8_t * cur  = &pv->mask_filtered[k][y*ref];
            uint8_t * dst = &pv->mask_temp[k][y*ref];
            
            int x;
            for( x = 1; x < width_limit; x++ )
            {
                dst[0] = 255;
                if( cur[0] == 255 )
                {
                    cur++;
                    dst++;
                    continue;
                }
                    
                count = 0;
                if( cur[-1-ref] == 255 )
                    ++count;
                if( cur[-ref] == 255 )
                    ++count;
                if( cur[+1-ref] == 255 )
                    ++count;
                if( cur[-1] == 255 )
                    ++count;
                if( cur[+1] == 255 )
                    ++count;
                if( cur[-1+ref] == 255 )
                    ++count;
                if( cur[+ref] == 255 )
                    ++count;
                if( cur[+1+ref] == 255 )
                    ++count;

                if( count < dilation_threshold )
                {
                    dst[0] = 0;
                }
                
                cur++;
                dst++;
            }
        }
    }
}


void filter_combing_mask( hb_filter_private_t * pv )
{
    int x, y, k;
    
    uint8_t * curp;
    uint8_t * cur;
    uint8_t * curn;
    uint8_t * dst;
    
    int h_count, v_count;
    
    int ref;
    
    int classic_filter = pv->filter_mode == FILTER_CLASSIC;
    
    int cur_is_combed = 0;
    
    for( k = 0; k < 1; k++ )
    {
        ref = pv->ref_stride[k];
        int width_limit = pv->width[k] - 1;
        int height_limit = pv->height[k] - 1;
        
        for( y = 0; y < height_limit; y++ )
        {
            curp = &pv->mask[k][(y-1)*ref];
            cur = &pv->mask[k][y*ref];
            curn = &pv->mask[k][(y+1)*ref];
            dst = ( classic_filter ) ? &pv->mask_filtered[k][y*ref] : &pv->mask_temp[k][y*ref] ;
            
            
            for( x = 0; x < width_limit; x++ )
            {
                h_count = 0;
                v_count = 0;
                
                if( cur[0] == 255 )
                {
                    if( x == 0 && cur[1] == 255 )
                    {
                        h_count++;
                    }
                    else if( x == width_limit && cur[-1] == 255 )
                    {
                        h_count++;
                    }
                    else if( cur[-1] == 255 && cur[1] == 255 )
                    {
                        h_count++;
                    }
                    
                    if( y == 0 && curn[0] == 255 )
                    {
                        v_count++;
                    }
                    else if( y == height_limit && curp[0] == 255 )
                    {
                        v_count++;
                    }
                    else if( curp[0] == 255 && curn[0] == 255 )
                    {
                        v_count++;
                    }
                    
                }
                
                if( classic_filter )
                {
                    if( h_count )
                        dst[0] = 255;
                    else
                        dst[0] = 0;
                }
                else
                {
                    if( h_count && v_count )
                        dst[0] = 255;
                    else
                        dst[0] = 0;
                }
                    
                curp++;
                cur++;
                curn++;
                dst++;
            }
        }
    }
}
*/


void erode_combing_mask( hb_filter_private_t * pv )
{
    //Take in mask, output to mask_temp
    int x, y, k;
    
    uint8_t * cur;
    uint8_t * dst;
    
    int count;
    
    int erosion_threshold = 2;
    
    int ref;
    for( k = 0; k < 1; k++ )
    {
        ref = pv->ref_stride[k];
        
        for( y = 1; y < pv->height[k] -1; y++ )
        {
            cur = &pv->mask_temp[k][y*ref];
            dst = &pv->mask_filtered[k][y*ref];
            
            for( x = 1; x < pv->width[k]-1; x++ )
            {
                if( cur[0] == 0 )
                {
                    dst[0] = 0;
                    cur++;
                    dst++;
                    continue;
                }
                    
                count = 0;
                if( cur[-ref-1] == 255 )
                    ++count;
                if( cur[-ref] == 255 )
                    ++count;
                if( cur[-ref+1] == 255 )
                    ++count;
                if( cur[-1] == 255 )
                    ++count;
                if( cur[+1] == 255 )
                    ++count;
                if( cur[ref-1] == 255 )
                    ++count;
                if( cur[ref] == 255 )
                    ++count;
                if( cur[ref+1] == 255 )
                    ++count;

                if( count < erosion_threshold )
                    dst[0] = 0;
                else
                    dst[0] = 255;
                    
                cur++;
                dst++;
            }
        }
    }
}

void dilate_combing_mask( hb_filter_private_t * pv )
{
    //Take in mask_temp, output to mask
    int x, y, k;
    
    uint8_t * curp;
    uint8_t * cur;
    uint8_t * curn;
    uint8_t * dst;
    
    int count;
    
    int dilation_threshold = 4;
    
    int ref;
    for( k = 0; k < 1; k++ )
    {
        ref = pv->ref_stride[k];
        
        for( y = 1; y < pv->height[k] -1; y++ )
        {
            curp = &pv->mask_filtered[k][(y-1)*ref];
            cur  = &pv->mask_filtered[k][y*ref];
            curn = &pv->mask_filtered[k][(y+1)*ref];
            dst = &pv->mask_temp[k][y*ref];
            
            for( x = 1; x < pv->width[k]-1; x++ )
            {
                if( cur[0] == 255 )
                {
                    dst[0] = 255;
                    curp++;
                    cur++;
                    curn++;
                    dst++;
                    continue;
                }
                    
                count = 0;
                if( curp[-1] == 255 )
                    ++count;
                if( curp[0] == 255 )
                    ++count;
                if( curp[+1] == 255 )
                    ++count;
                if( cur[-1] == 255 )
                    ++count;
                if( cur[+1] == 255 )
                    ++count;
                if( curn[-1] == 255 )
                    ++count;
                if( curn[0] == 255 )
                    ++count;
                if( curn[+1] == 255 )
                    ++count;

                if( count >= dilation_threshold )
                    dst[0] = 255;
                else
                    dst[0] = 0;
                    
                curp++;
                cur++;
                curn++;
                dst++;
            }
        }
    }
}

void filter_combing_mask( hb_filter_private_t * pv )
{
    int x, y, k;
    
    uint8_t * curp;
    uint8_t * cur;
    uint8_t * curn;
    uint8_t * dst;
    
    int h_count, v_count;
    
    int ref;
    for( k = 0; k < 1; k++ )
    {
        ref = pv->ref_stride[k];
        
        for( y = 0; y < pv->height[k] -1; y++ )
        {
            curp = &pv->mask[k][(y-1)*ref];
            cur = &pv->mask[k][y*ref];
            curn = &pv->mask[k][(y+1)*ref];
            dst = (pv->filter_mode == FILTER_CLASSIC ) ? &pv->mask_filtered[k][y*ref] : &pv->mask_temp[k][y*ref] ;
            
            for( x = 0; x < pv->width[k]-1; x++ )
            {
                
                h_count = v_count = 0;
                if( x == 0 )
                {
                    if( cur[0] == 255 && cur[1] == 255 )
                        h_count++;
                }
                else if( x == pv->width[k]-1 )
                {
                    if( cur[-1] == 255 && cur[0] == 255 )
                        h_count++;
                }
                else
                {
                    if(cur[-1] == 255 && cur[0] == 255 && cur[1] == 255 )
                        h_count++;
                }
                
                if( y == 0 )
                {
                    if( cur[0] == 255 && curn[0] == 255 )
                        v_count++;
                }
                else if( y == pv->height[k]-1 )
                {
                    if( curp[0] == 255 && cur[0] == 255 )
                        v_count++;
                }
                else
                {
                    if( curp[0] == 255 && cur[0] == 255 && curn[0] == 255 )
                        v_count++;
                }
                
                if( h_count & pv->filter_mode == FILTER_CLASSIC )
                    dst[0] = 255;
                else if( pv->filter_mode == FILTER_CLASSIC )
                    dst[0] = 0;
                else if( h_count && v_count )
                    dst[0] = 255;
                else
                    dst[0] = 0;
                    
                curp++;
                cur++;
                curn++;
                dst++;
            }
        }
    }
}

int check_filtered_combing_mask( hb_filter_private_t * pv )
{
    /* Go through the mask in X*Y blocks. If any of these windows
       have threshold or more combed pixels, consider the whole
       frame to be combed and send it on to be deinterlaced.     */

    /* Block mask threshold -- The number of pixels
       in a block_width * block_height window of
       he mask that need to show combing for the
       whole frame to be seen as such.            */
    int threshold       = pv->block_threshold;
    int block_width     = pv->block_width;
    int block_height    = pv->block_height;
    int block_x, block_y;
    int block_score = 0; int send_to_blend = 0;
    uint8_t * mask_p;
    int x, y, k;

    for( k = 0; k < 1; k++ )
    {
        int ref_stride = pv->ref_stride[k];
        pv->mask_box_x = -1;
        pv->mask_box_y = -1;
        pv->mask_box_color = 0;
        
        for( y = 0; y < ( pv->height[k] - block_height ); y = y + block_height )
        {
            for( x = 0; x < ( pv->width[k] - block_width ); x = x + block_width )
            {
                block_score = 0;
                
                for( block_y = 0; block_y < block_height; block_y++ )
                {
                    int mask_y = y + block_y;
                    mask_p = &pv->mask_filtered[k][mask_y*ref_stride + x];
                    
                    for( block_x = 0; block_x < block_width; block_x++ )
                    {

                        if( mask_p[ 0 ] == 255 )
                            block_score++;
                        mask_p++;
                    }
                }

                if( block_score >= ( threshold / 2 ) )
                {
#if 0
                    hb_log("decomb: frame %i | score %i | type %s", pv->deinterlaced_frames + pv->blended_frames +  pv->unfiltered_frames + 1, block_score, pv->buf_settings->flags & 16 ? "Film" : "Video");
#endif
                    pv->mask_box_x = x;
                    pv->mask_box_y = y;

                   if ( block_score <= threshold && !( pv->buf_settings->flags & 16) )
                    {
                        /* Blend video content that scores between
                           ( threshold / 2 ) and threshold.        */
                        send_to_blend = 1;
                        pv->mask_box_color = 2;
                    }
                    else if( block_score > threshold )
                    {
                        if( pv->buf_settings->flags & 16 )
                        {
                            /* Blend progressive content above the threshold.*/
                            pv->mask_box_color = 2;
                            return 2;
                        }
                        else
                        {
                            /* Yadif deinterlace video content above the threshold. */
                            pv->mask_box_color = 1;
                            return 1;
                        }
                    }
                }
            }
        } 
    }
    
    if( send_to_blend )
    {
        return 2;
    }
    else
    {
        /* Consider this frame to be uncombed. */
        return 0;
    }
}

int check_combing_mask( hb_filter_private_t * pv )
{
    /* Go through the mask in X*Y blocks. If any of these windows
       have threshold or more combed pixels, consider the whole
       frame to be combed and send it on to be deinterlaced.     */

    /* Block mask threshold -- The number of pixels
       in a block_width * block_height window of
       he mask that need to show combing for the
       whole frame to be seen as such.            */
    int threshold       = pv->block_threshold;
    int block_width     = pv->block_width;
    int block_height    = pv->block_height;
    int block_x, block_y;
    int block_score = 0; int send_to_blend = 0;
    uint8_t * mask_p;
    int x, y, k;

    for( k = 0; k < 1; k++ )
    {
        int ref_stride = pv->ref_stride[k];
        for( y = 0; y < ( pv->height[k] - block_height ); y = y + block_height )
        {
            for( x = 0; x < ( pv->width[k] - block_width ); x = x + block_width )
            {
                block_score = 0;
                
                for( block_y = 0; block_y < block_height; block_y++ )
                {
                    int mask_y = y + block_y;
                    mask_p = &pv->mask[k][mask_y*ref_stride + x];
                    
                    for( block_x = 0; block_x < block_width; block_x++ )
                    {
                        /* We only want to mark a pixel in a block as combed
                           if the adjacent pixels are as well. Got to
                           handle the sides separately.       */
                        if( (x + block_x) == 0 )
                        {
                            if( mask_p[ 0 ] == 255 &&
                                mask_p[ 1 ] == 255 )
                                    block_score++;
                        }
                        else if( (x + block_x) == (pv->width[k] -1) )
                        {
                            if( mask_p[ -1 ] == 255 &&
                                mask_p[  0 ] == 255 )
                                    block_score++;
                        }
                        else
                        {
                            if( mask_p[ -1 ] == 255 &&
                                mask_p[  0 ] == 255 &&
                                mask_p[  1 ] == 255 )
                                    block_score++;
                        }

//                           if( (y + block_y) == 0 )
//                           {
//                               if( mask_p[ 0 ] == 255 &&
//                                   mask_p[ ref_stride ] == 255 )
//                                       block_score++;
//                           }
//                           else if( (y + block_y) == (pv->height[k] -1) )
//                           {
//                               if( mask_p[ -ref_stride ] == 255 &&
//                                   mask_p[  0 ] == 255 )
//                                       block_score++;
//                           }
//                           else
//                           {
//                               if( mask_p[ -ref_stride ] == 255 &&
//                                   mask_p[  0 ] == 255 &&
//                                   mask_p[  +ref_stride ] == 255 )
//                                       block_score++;
//                           }
//
                        mask_p++;
                    }
                }

                if( block_score >= ( threshold / 2 ) )
                {
#if 0
                    hb_log("decomb: frame %i | score %i | type %s", pv->deinterlaced_frames + pv->blended_frames +  pv->unfiltered_frames + 1, block_score, pv->buf_settings->flags & 16 ? "Film" : "Video");
#endif

                    pv->mask_box_x = x;
                    pv->mask_box_y = y;

                    if ( block_score <= threshold && !( pv->buf_settings->flags & 16) )
                    {
                        /* Blend video content that scores between
                           ( threshold / 2 ) and threshold.        */
                        send_to_blend = 1;
                        pv->mask_box_color = 2;
                    }
                    else if( block_score > threshold )
                    {
                        if( pv->buf_settings->flags & 16 )
                        {
                            /* Blend progressive content above the threshold.*/
                            pv->mask_box_color = 2;
                            return 2;
                        }
                        else
                        {
                            /* Yadif deinterlace video content above the threshold. */
                            pv->mask_box_color = 1;
                            return 1;
                        }
                    }
                }
            }
        } 
    }
    
    if( send_to_blend )
    {
        return 2;
    }
    else
    {
        /* Consider this frame to be uncombed. */
        return 0;
    }
}

void build_gamma_lut( hb_filter_private_t * pv )
{
    int i;
    for( i = 0; i < 256; i++ )
    {
        pv->gamma_lut[i] = pow( ( (float)i / (float)255 ), 2.2f );
    }
}

float scale_gamma( int pixel, hb_filter_private_t * pv )
{
    return pv->gamma_lut[pixel];
}

#if 0
void detect_gamma_combed_segment( hb_filter_private_t * pv, int segment_start, int segment_stop )
{
    /* A mish-mash of various comb detection tricks
       picked up from neuron2's Decomb plugin for
       AviSynth and tritical's IsCombedT and
       IsCombedTIVTC plugins.                       */
       
    int x, y, k, width, height;
    int first_frame = ( pv->deinterlaced_frames==0 && pv->blended_frames==0 && pv->unfiltered_frames==0);
    /* Comb scoring algorithm */
//    int spatial_metric  = pv->spatial_metric;
    /* Motion threshold */
    const float mthresh         = (float)pv->motion_threshold / (float)255;
    /* Spatial threshold */
    const float athresh         = (float)pv->spatial_threshold / (float)255;
//    float athresh_squared = athresh * athresh;
    const float athresh6        = 6 *athresh;

    static const float __attribute__ ((aligned (16))) gamma_lut[256] = { 0.000000,0.000005,0.000023,0.000057,0.000107,0.000175,0.000262,0.000367,0.000493,0.000638,0.000805,0.000992,0.001202,0.001433,0.001687,0.001963,0.002263,0.002586,0.002932,0.003303,0.003697,0.004116,0.004560,0.005028,0.005522,0.006041,0.006585,0.007155,0.007751,0.008373,0.009021,0.009696,0.010398,0.011126,0.011881,0.012664,0.013473,0.014311,0.015175,0.016068,0.016988,0.017936,0.018913,0.019918,0.020951,0.022013,0.023104,0.024223,0.025371,0.026549,0.027755,0.028991,0.030257,0.031551,0.032876,0.034230,0.035614,0.037029,0.038473,0.039947,0.041452,0.042987,0.044553,0.046149,0.047776,0.049433,0.051122,0.052842,0.054592,0.056374,0.058187,0.060032,0.061907,0.063815,0.065754,0.067725,0.069727,0.071761,0.073828,0.075926,0.078057,0.080219,0.082414,0.084642,0.086901,0.089194,0.091518,0.093876,0.096266,0.098689,0.101145,0.103634,0.106156,0.108711,0.111299,0.113921,0.116576,0.119264,0.121986,0.124741,0.127530,0.130352,0.133209,0.136099,0.139022,0.141980,0.144972,0.147998,0.151058,0.154152,0.157281,0.160444,0.163641,0.166872,0.170138,0.173439,0.176774,0.180144,0.183549,0.186989,0.190463,0.193972,0.197516,0.201096,0.204710,0.208360,0.212044,0.215764,0.219520,0.223310,0.227137,0.230998,0.234895,0.238828,0.242796,0.246800,0.250840,0.254916,0.259027,0.263175,0.267358,0.271577,0.275833,0.280124,0.284452,0.288816,0.293216,0.297653,0.302126,0.306635,0.311181,0.315763,0.320382,0.325037,0.329729,0.334458,0.339223,0.344026,0.348865,0.353741,0.358654,0.363604,0.368591,0.373615,0.378676,0.383775,0.388910,0.394083,0.399293,0.404541,0.409826,0.415148,0.420508,0.425905,0.431340,0.436813,0.442323,0.447871,0.453456,0.459080,0.464741,0.470440,0.476177,0.481952,0.487765,0.493616,0.499505,0.505432,0.511398,0.517401,0.523443,0.529523,0.535642,0.541798,0.547994,0.554227,0.560499,0.566810,0.573159,0.579547,0.585973,0.592438,0.598942,0.605484,0.612066,0.618686,0.625345,0.632043,0.638779,0.645555,0.652370,0.659224,0.666117,0.673049,0.680020,0.687031,0.694081,0.701169,0.708298,0.715465,0.722672,0.729919,0.737205,0.744530,0.751895,0.759300,0.766744,0.774227,0.781751,0.789314,0.796917,0.804559,0.812241,0.819964,0.827726,0.835528,0.843370,0.851252,0.859174,0.867136,0.875138,0.883180,0.891262,0.899384,0.907547,0.915750,0.923993,0.932277,0.940601,0.948965,0.957370,0.965815,0.974300,0.982826,0.991393,1.000000 };

    /* One pas for Y, one pass for U, one pass for V */    
    for( k = 0; k < 1; k++ )
    {
        int ref_stride  = pv->ref_stride[k];
        width           = pv->width[k];
        height          = pv->height[k];
        
        /* Comb detection has to start at y = 2 and end at
           y = height - 2, because it needs to examine
           2 pixels above and 2 below the current pixel.      */
        if( segment_start < 2 )
            segment_start = 2;
        if( segment_stop > height - 2 )
            segment_stop = height - 2;

        for( y =  segment_start; y < segment_stop; y++ )
        {
            /* These are just to make the buffer locations easier to read. */
            int up_2    = -2*ref_stride ;
            int up_1    = -1*ref_stride;
            int down_1 = ref_stride;
            int down_2 = 2*ref_stride;
            
            /* We need to examine a column of 5 pixels
               in the prev, cur, and next frames.      */
            uint8_t * cur = &pv->ref[1][k][y*ref_stride];
            uint8_t * prev = &pv->ref[0][k][y*ref_stride];
            uint8_t * next = &pv->ref[2][k][y*ref_stride];
            uint8_t * mask = &pv->mask[k][y*ref_stride];
            
            for( x = 0; x < width; x++ )
            {
                const float up_diff = gamma_lut[cur[0]] - gamma_lut[cur[-ref_stride]];
                const float down_diff = gamma_lut[cur[0]] - gamma_lut[cur[+ref_stride]];
                
                mask[0] = 0;
                if( ( up_diff >  athresh && down_diff >  athresh ) ||
                    ( up_diff < -athresh && down_diff < -athresh ) )
                {
                    /* The pixel above and below are different,
                       and they change in the same "direction" too.*/
                    int motion = 1; // Assume user does not want motion detect
                    if( mthresh > 0 )
                    {
                        motion = 0; // User does want motion detect
                        /* Make sure there's sufficient motion between frame t-1 to frame t+1. */
                        if( fabs( gamma_lut[prev[0]]      - gamma_lut[cur[0]] ) > mthresh &&
                            fabs( gamma_lut[cur[up_1]]    - gamma_lut[next[up_1]]    ) > mthresh &&
                            fabs( gamma_lut[cur[down_1]]  - gamma_lut[next[down_1]]    ) > mthresh )
                                motion++;
                        if( fabs( gamma_lut[next[0]]      - gamma_lut[cur[0]] ) > mthresh &&
                            fabs( gamma_lut[prev[up_1]]   - gamma_lut[cur[up_1]] ) > mthresh &&
                            fabs( gamma_lut[prev[down_1]] - gamma_lut[cur[down_1]] ) > mthresh )
                                motion++;
                                
//                            hb_log("prev->cur motion: %f, mthresh: %f", fabs( scale_gamma( prev[0] ) - scale_gamma( cur[0] ) ), mthresh);
                    }
                           
                    
                    if( motion || first_frame )
                    {

                        /* Tritical's noise-resistant combing scorer.
                           The check is done on a bob+blur convolution. */
                        float combing = fabs( gamma_lut[cur[up_2]]
                                         + ( 4 * gamma_lut[cur[0]] )
                                         + gamma_lut[cur[down_2]]
                                         - ( 3 * ( gamma_lut[cur[up_1]]
                                                 + gamma_lut[cur[down_1]] ) ) );
                        /* If the frame is sufficiently combed,
                           then mark it down on the mask as 255. */
                           
 //                              hb_log("combing: %f, athresh6: %f", combing, athresh6);
                        if( combing > athresh6 )
                        {
                            mask[0] = 255;
                        }
                    }
                }
                
                cur++;
                prev++;
                next++;
                mask++;
            }
        }
    }
}
#endif

void detect_gamma_combed_segment( hb_filter_private_t * pv, int segment_start, int segment_stop )
{
    /* A mish-mash of various comb detection tricks
       picked up from neuron2's Decomb plugin for
       AviSynth and tritical's IsCombedT and
       IsCombedTIVTC plugins.                       */
       
    int x, y, k, width, height;
    
    /* Comb scoring algorithm */
    float spatial_metric  = (float)pv->spatial_metric / (float)255;
    /* Motion threshold */
    float mthresh         = (float)pv->motion_threshold / (float)255;
    /* Spatial threshold */
    float athresh         = (float)pv->spatial_threshold / (float)255;
    float athresh_squared = athresh * athresh;
    float athresh6        = 6 *athresh;

    /* One pas for Y, one pass for U, one pass for V */    
    for( k = 0; k < 1; k++ )
    {
        int ref_stride  = pv->ref_stride[k];
        width           = pv->width[k];
        height          = pv->height[k];
        
        /* Comb detection has to start at y = 2 and end at
           y = height - 2, because it needs to examine
           2 pixels above and 2 below the current pixel.      */
        if( segment_start < 2 )
            segment_start = 2;
        if( segment_stop > height - 2 )
            segment_stop = height - 2;
            
        for( y =  segment_start; y < segment_stop; y++ )
        {
            /* These are just to make the buffer locations easier to read. */
            int up_2    = -2*ref_stride ;
            int up_1    = -1*ref_stride;
            int down_1 = ref_stride;
            int down_2 = 2*ref_stride;
            
            /* We need to examine a column of 5 pixels
               in the prev, cur, and next frames.      */
            uint8_t * cur = &pv->ref[1][k][y*ref_stride];
            uint8_t * prev = &pv->ref[0][k][y*ref_stride];
            uint8_t * next = &pv->ref[2][k][y*ref_stride];
            uint8_t * mask = &pv->mask[k][y*ref_stride];
            
            for( x = 0; x < width; x++ )
            {
                float up_diff   = pv->gamma_lut[cur[0]] - pv->gamma_lut[cur[up_1]];
                float down_diff = pv->gamma_lut[cur[0]] - pv->gamma_lut[cur[down_1]];
                
                if( ( up_diff >  athresh && down_diff >  athresh ) ||
                    ( up_diff < -athresh && down_diff < -athresh ) )
                {
                    /* The pixel above and below are different,
                       and they change in the same "direction" too.*/
                    int motion = 0;
                    if( mthresh > 0 )
                    {
                        /* Make sure there's sufficient motion between frame t-1 to frame t+1. */
                        if( fabs( pv->gamma_lut[prev[0]]      - pv->gamma_lut[cur[0]] ) > mthresh &&
                            fabs( pv->gamma_lut[cur[up_1]]    - pv->gamma_lut[next[up_1]]    ) > mthresh &&
                            fabs( pv->gamma_lut[cur[down_1]]  - pv->gamma_lut[next[down_1]]    ) > mthresh )
                                motion++;
                        if( fabs( pv->gamma_lut[next[0]]      - pv->gamma_lut[cur[0]] ) > mthresh &&
                            fabs( pv->gamma_lut[prev[up_1]]   - pv->gamma_lut[cur[up_1]] ) > mthresh &&
                            fabs( pv->gamma_lut[prev[down_1]] - pv->gamma_lut[cur[down_1]] ) > mthresh )
                                motion++;
                                
//                            hb_log("prev->cur motion: %f, mthresh: %f", fabs( scale_gamma( prev[0] ) - scale_gamma( cur[0] ) ), mthresh);
                    }
                    else
                    {
                        /* User doesn't want to check for motion,
                           so move on to the spatial check.       */
                        motion = 1;
                    }
                           
                    if( motion || ( pv->deinterlaced_frames==0 && pv->blended_frames==0 && pv->unfiltered_frames==0) )
                    {

                        /* Tritical's noise-resistant combing scorer.
                           The check is done on a bob+blur convolution. */
                        float combing = fabs( pv->gamma_lut[cur[up_2]]
                                         + ( 4 * pv->gamma_lut[cur[0]] )
                                         + pv->gamma_lut[cur[down_2]]
                                         - ( 3 * ( pv->gamma_lut[cur[up_1]]
                                                 + pv->gamma_lut[cur[down_1]] ) ) );
                        /* If the frame is sufficiently combed,
                           then mark it down on the mask as 255. */
                           
 //                              hb_log("combing: %f, athresh6: %f", combing, athresh6);
                        if( combing > athresh6 )
                        {
                            mask[0] = 255;
                        }
                        else
                        {
                            mask[0] = 0;
                        }

                    }
                    else
                    {
                        mask[0] = 0;
                    }
                }
                else
                {
                    mask[0] = 0;
                }
                
                cur++;
                prev++;
                next++;
                mask++;
            }
        }
    }
}


void detect_combed_segment( hb_filter_private_t * pv, int segment_start, int segment_stop )
{
    /* A mish-mash of various comb detection tricks
       picked up from neuron2's Decomb plugin for
       AviSynth and tritical's IsCombedT and
       IsCombedTIVTC plugins.                       */
       
    int x, y, k, width, height;
    
    /* Comb scoring algorithm */
    int spatial_metric  = pv->spatial_metric;
    /* Motion threshold */
    int mthresh         = pv->motion_threshold;
    /* Spatial threshold */
    int athresh         = pv->spatial_threshold;
    int athresh_squared = athresh * athresh;
    int athresh6        = 6 *athresh;

    /* One pas for Y, one pass for U, one pass for V */    
    for( k = 0; k < 1; k++ )
    {
        int ref_stride  = pv->ref_stride[k];
        width           = pv->width[k];
        height          = pv->height[k];
        
        /* Comb detection has to start at y = 2 and end at
           y = height - 2, because it needs to examine
           2 pixels above and 2 below the current pixel.      */
        if( segment_start < 2 )
            segment_start = 2;
        if( segment_stop > height - 2 )
            segment_stop = height - 2;
            
        for( y =  segment_start; y < segment_stop; y++ )
        {
            /* These are just to make the buffer locations easier to read. */
            int up_2    = -2*ref_stride ;
            int up_1    = -1*ref_stride;
            int down_1 = ref_stride;
            int down_2 = 2*ref_stride;
            
            /* We need to examine a column of 5 pixels
               in the prev, cur, and next frames.      */
            uint8_t * cur = &pv->ref[1][k][y*ref_stride];
            uint8_t * prev = &pv->ref[0][k][y*ref_stride];
            uint8_t * next = &pv->ref[2][k][y*ref_stride];
            uint8_t * mask = &pv->mask[k][y*ref_stride];
            
            for( x = 0; x < width; x++ )
            {
                int up_diff = cur[0] - cur[up_1];
                int down_diff = cur[0] - cur[down_1];
                
                if( ( up_diff >  athresh && down_diff >  athresh ) ||
                    ( up_diff < -athresh && down_diff < -athresh ) )
                {
                    /* The pixel above and below are different,
                       and they change in the same "direction" too.*/
                    int motion = 0;
                    if( mthresh > 0 )
                    {
                        /* Make sure there's sufficient motion between frame t-1 to frame t+1. */
                        if( abs( prev[0] - cur[0] ) > mthresh &&
                            abs(  cur[up_1] - next[up_1]    ) > mthresh &&
                            abs(  cur[down_1] - next[down_1]    ) > mthresh )
                                motion++;
                        if( abs(     next[0] - cur[0] ) > mthresh &&
                            abs( prev[up_1] - cur[up_1] ) > mthresh &&
                            abs( prev[down_1] - cur[down_1] ) > mthresh )
                                motion++;
                    }
                    else
                    {
                        /* User doesn't want to check for motion,
                           so move on to the spatial check.       */
                        motion = 1;
                    }
                           
                    if( motion || ( pv->deinterlaced_frames==0 && pv->blended_frames==0 && pv->unfiltered_frames==0) )
                    {
                           /* That means it's time for the spatial check.
                              We've got several options here.             */
                        if( spatial_metric == 0 )
                        {
                            /* Simple 32detect style comb detection */
                            if( ( abs( cur[0] - cur[down_2] ) < 10  ) &&
                                ( abs( cur[0] - cur[down_1] ) > 15 ) )
                            {
                                mask[0] = 255;
                            }
                            else
                            {
                                mask[0] = 0;
                            }
                        }
                        else if( spatial_metric == 1 )
                        {
                            /* This, for comparison, is what IsCombed uses.
                               It's better, but still noise senstive.      */
                               int combing = ( cur[up_1] - cur[0] ) *
                                             ( cur[down_1] - cur[0] );
                               
                               if( combing > athresh_squared )
                                   mask[0] = 255; 
                               else
                                   mask[0] = 0;
                        }
                        else if( spatial_metric == 2 )
                        {
                            /* Tritical's noise-resistant combing scorer.
                               The check is done on a bob+blur convolution. */
                            int combing = abs( cur[up_2]
                                             + ( 4 * cur[0] )
                                             + cur[down_2]
                                             - ( 3 * ( cur[up_1]
                                                     + cur[down_1] ) ) );

                            /* If the frame is sufficiently combed,
                               then mark it down on the mask as 255. */
                            if( combing > athresh6 )
                            {
                                mask[0] = 255;
                            }
                            else
                            {
                                mask[0] = 0;
                            }
                        }
                    }
                    else
                    {
                        mask[0] = 0;
                    }
                }
                else
                {
                    mask[0] = 0;
                }
                
                cur++;
                prev++;
                next++;
                mask++;
            }
        }
    }
}

// This function calls all the eedi2 filters in sequence for a given plane.
// It outputs the final interpolated image to pv->eedi_full[DST2PF].
void eedi2_interpolate_plane( hb_filter_private_t * pv, int k )
{
    /* We need all these pointers. No, seriously.
       I swear. It's not a joke. They're used.
       All nine of them.                         */
    uint8_t * mskp = pv->eedi_half[MSKPF][k];
    uint8_t * srcp = pv->eedi_half[SRCPF][k];
    uint8_t * tmpp = pv->eedi_half[TMPPF][k];
    uint8_t * dstp = pv->eedi_half[DSTPF][k];
    uint8_t * dst2p = pv->eedi_full[DST2PF][k];
    uint8_t * tmp2p2 = pv->eedi_full[TMP2PF2][k];
    uint8_t * msk2p = pv->eedi_full[MSK2PF][k];
    uint8_t * tmp2p = pv->eedi_full[TMP2PF][k];
    uint8_t * dst2mp = pv->eedi_full[DST2MPF][k];
    int * cx2 = pv->cx2;
    int * cy2 = pv->cy2;
    int * cxy = pv->cxy;
    int * tmpc = pv->tmpc;

    int pitch = pv->ref_stride[k];
    int height = pv->height[k]; int width = pv->width[k];
    int half_height = height / 2;

    // edge mask
    eedi2_build_edge_mask( mskp, pitch, srcp, pitch,
                     pv->magnitude_threshold, pv->variance_threshold, pv->laplacian_threshold, 
                     half_height, width );
    eedi2_erode_edge_mask( mskp, pitch, tmpp, pitch, pv->erosion_threshold, half_height, width );
    eedi2_dilate_edge_mask( tmpp, pitch, mskp, pitch, pv->dilation_threshold, half_height, width );
    eedi2_erode_edge_mask( mskp, pitch, tmpp, pitch, pv->erosion_threshold, half_height, width );
    eedi2_remove_small_gaps( tmpp, pitch, mskp, pitch, half_height, width );

    // direction mask
    eedi2_calc_directions( k, mskp, pitch, srcp, pitch, tmpp, pitch,
                     pv->maximum_search_distance, pv->noise_threshold,
                     half_height, width );
    eedi2_filter_dir_map( mskp, pitch, tmpp, pitch, dstp, pitch, half_height, width );
    eedi2_expand_dir_map( mskp, pitch, dstp, pitch, tmpp, pitch, half_height, width );
    eedi2_filter_map( mskp, pitch, tmpp, pitch, dstp, pitch, half_height, width );

    // upscale 2x vertically
    eedi2_upscale_by_2( srcp, dst2p, half_height, pitch );
    eedi2_upscale_by_2( dstp, tmp2p2, half_height, pitch );
    eedi2_upscale_by_2( mskp, msk2p, half_height, pitch );

    // upscale the direction mask
    eedi2_mark_directions_2x( msk2p, pitch, tmp2p2, pitch, tmp2p, pitch, pv->tff, height, width );
    eedi2_filter_dir_map_2x( msk2p, pitch, tmp2p, pitch,  dst2mp, pitch, pv->tff, height, width );
    eedi2_expand_dir_map_2x( msk2p, pitch, dst2mp, pitch, tmp2p, pitch, pv->tff, height, width );
    eedi2_fill_gaps_2x( msk2p, pitch, tmp2p, pitch, dst2mp, pitch, pv->tff, height, width );
    eedi2_fill_gaps_2x( msk2p, pitch, dst2mp, pitch, tmp2p, pitch, pv->tff, height, width );

    // interpolate a full-size plane
    eedi2_interpolate_lattice( k, tmp2p, pitch, dst2p, pitch, tmp2p2, pitch, pv->tff,
                         pv->noise_threshold, height, width );

    if( pv->post_processing == 1 || pv->post_processing == 3 )
    {
        // make sure the edge directions are consistent
        eedi2_bit_blit( tmp2p2, pitch, tmp2p, pitch, pv->width[k], pv->height[k] );
        eedi2_filter_dir_map_2x( msk2p, pitch, tmp2p, pitch, dst2mp, pitch, pv->tff, height, width );
        eedi2_expand_dir_map_2x( msk2p, pitch, dst2mp, pitch, tmp2p, pitch, pv->tff, height, width );
        eedi2_post_process( tmp2p, pitch, tmp2p2, pitch, dst2p, pitch, pv->tff, height, width );
    }
    if( pv->post_processing == 2 || pv->post_processing == 3 )
    {
        // filter junctions and corners
        eedi2_gaussian_blur1( srcp, pitch, tmpp, pitch, srcp, pitch, half_height, width );
        eedi2_calc_derivatives( srcp, pitch, half_height, width, cx2, cy2, cxy );
        eedi2_gaussian_blur_sqrt2( cx2, tmpc, cx2, pitch, half_height, width);
        eedi2_gaussian_blur_sqrt2( cy2, tmpc, cy2, pitch, half_height, width);
        eedi2_gaussian_blur_sqrt2( cxy, tmpc, cxy, pitch, half_height, width);
        eedi2_post_process_corner( cx2, cy2, cxy, pitch, tmp2p2, pitch, dst2p, pitch, height, width, pv->tff );
    }
}

/*
 *  eedi2 interpolate this plane in a single thread.
 */
void eedi2_filter_thread( void *thread_args_v )
{
    eedi2_arguments_t *eedi2_work = NULL;
    hb_filter_private_t * pv;
    int run = 1;
    int plane;
    eedi2_thread_arg_t *thread_args = thread_args_v;

    pv = thread_args->pv;
    plane = thread_args->plane;

    hb_log("eedi2 thread started for plane %d", plane);

    while( run )
    {
        /*
         * Wait here until there is work to do. hb_lock() blocks until
         * render releases it to say that there is more work to do.
         */
        hb_lock( pv->eedi2_begin_lock[plane] );

        eedi2_work = &pv->eedi2_arguments[plane];

        if( eedi2_work->stop )
        {
            /*
             * No more work to do, exit this thread.
             */
            run = 0;
            continue;
        } 

        /*
         * Process plane
         */
            eedi2_interpolate_plane( pv, plane );
        
        /*
         * Finished this segment, let everyone know.
         */
        hb_unlock( pv->eedi2_complete_lock[plane] );
    }
    free( thread_args_v );
}

// Sets up the input field planes for EEDI2 in pv->eedi_half[SRCPF]
// and then runs eedi2_filter_thread for each plane.
void eedi2_planer( hb_filter_private_t * pv )
{
    /* Copy the first field from the source to a half-height frame. */
    int i;
    for( i = 0;  i < 3; i++ )
    {
        int pitch = pv->ref_stride[i];
        int start_line = !pv->tff;
        eedi2_fill_half_height_buffer_plane( &pv->ref[1][i][pitch*start_line], pv->eedi_half[SRCPF][i], pitch, pv->height[i] );
    }
    
    int plane;
    for( plane = 0; plane < 3; plane++ )
    {  
        /*
         * Let the thread for this plane know that we've setup work 
         * for it by releasing the begin lock (ensuring that the
         * complete lock is already locked so that we block when
         * we try to lock it again below).
         */
        hb_lock( pv->eedi2_complete_lock[plane] );
        hb_unlock( pv->eedi2_begin_lock[plane] );
    }

    /*
     * Wait until all three threads have completed by trying to get
     * the complete lock that we locked earlier for each thread, which
     * will block until that thread has completed the work on that
     * plane.
     */
    for( plane = 0; plane < 3; plane++ )
    {
        hb_lock( pv->eedi2_complete_lock[plane] );
        hb_unlock( pv->eedi2_complete_lock[plane] );
    }
}


/*
 * comb detect this segment of all three planes in a single thread.
 */
void decomb_filter_thread( void *thread_args_v )
{
    decomb_arguments_t *decomb_work = NULL;
    hb_filter_private_t * pv;
    int run = 1;
    int segment, segment_start, segment_stop, plane;
    decomb_thread_arg_t *thread_args = thread_args_v;

    pv = thread_args->pv;
    segment = thread_args->segment;

    hb_log("decomb thread started for segment %d", segment);

    while( run )
    {
        /*
         * Wait here until there is work to do. hb_lock() blocks until
         * render releases it to say that there is more work to do.
         */
        hb_lock( pv->decomb_begin_lock[segment] );

        decomb_work = &pv->decomb_arguments[segment];

        if( decomb_work->stop )
        {
            /*
             * No more work to do, exit this thread.
             */
            run = 0;
            continue;
        } 

        /*
         * Process segment (for now just from luma)
         */
        for( plane = 0; plane < 1; plane++)
        {

            int h = pv->height[plane];
            segment_start = ( h / pv->cpu_count ) * segment;
            if( segment == pv->cpu_count - 1 )
            {
                /*
                 * Final segment
                 */
                segment_stop = h;
            } else {
                segment_stop = ( h / pv->cpu_count ) * ( segment + 1 );
            }
            
            if( pv->mode & MODE_GAMMA )
            {
                detect_gamma_combed_segment( pv, segment_start, segment_stop );
            }
            else
            {
                detect_combed_segment( pv, segment_start, segment_stop );
            }
        }
        /*
         * Finished this segment, let everyone know.
         */
        hb_unlock( pv->decomb_complete_lock[segment] );
    }
    free( thread_args_v );
}

int comb_segmenter( hb_filter_private_t * pv )
{
    int segment;

    for( segment = 0; segment < pv->cpu_count; segment++ )
    {  
        /*
         * Let the thread for this plane know that we've setup work 
         * for it by releasing the begin lock (ensuring that the
         * complete lock is already locked so that we block when
         * we try to lock it again below).
         */
        hb_lock( pv->decomb_complete_lock[segment] );
        hb_unlock( pv->decomb_begin_lock[segment] );
    }

    /*
     * Wait until all three threads have completed by trying to get
     * the complete lock that we locked earlier for each thread, which
     * will block until that thread has completed the work on that
     * plane.
     */
    for( segment = 0; segment < pv->cpu_count; segment++ )
    {
        hb_lock( pv->decomb_complete_lock[segment] );
        hb_unlock( pv->decomb_complete_lock[segment] );
    }
    
    if( pv->mode & MODE_FILTER )
    {
        filter_combing_mask( pv );
        if( pv->filter_mode == FILTER_ERODE_DILATE )
        {
            erode_combing_mask( pv );
            dilate_combing_mask( pv );
            erode_combing_mask( pv );
        }
        return check_filtered_combing_mask( pv );
    }
    else
    {
        return check_combing_mask( pv );
    }
}

static void yadif_filter_line( uint8_t *dst,
                               uint8_t *prev,
                               uint8_t *cur,
                               uint8_t *next,
                               int plane,
                               int parity,
                               int y,
                               hb_filter_private_t * pv )
{
    /* While prev and next point to the previous and next frames,
       prev2 and next2 will shift depending on the parity, usually 1.
       They are the previous and next fields, the fields temporally adjacent
       to the other field in the current frame--the one not being filtered.  */
    uint8_t *prev2 = parity ? prev : cur ;
    uint8_t *next2 = parity ? cur  : next;
    
    int w = pv->width[plane];
    int refs = pv->ref_stride[plane];
    int x;
    int eedi2_mode = ( pv->mode & MODE_EEDI2 );
    
    /* We can replace spatial_pred with this interpolation*/
    uint8_t * eedi2_guess = &pv->eedi_full[DST2PF][plane][y*refs];

    /* Decomb's cubic interpolation can only function when there are
       three samples above and below, so regress to yadif's traditional
       two-tap interpolation when filtering at the top and bottom edges. */
    int vertical_edge = 0;
    if( ( y < 3 ) || ( y > ( pv->height[plane] - 4 ) )  )
        vertical_edge = 1;

    for( x = 0; x < w; x++)
    {
        /* Pixel above*/
        int c              = cur[-refs];
        /* Temporal average: the current location in the adjacent fields */
        int d              = (prev2[0] + next2[0])>>1;
        /* Pixel below */
        int e              = cur[+refs];
        
        /* How the current pixel changes between the adjacent fields */
        int temporal_diff0 = ABS(prev2[0] - next2[0]);
        /* The average of how much the pixels above and below change from the frame before to now. */
        int temporal_diff1 = ( ABS(prev[-refs] - cur[-refs]) + ABS(prev[+refs] - cur[+refs]) ) >> 1;
        /* The average of how much the pixels above and below change from now to the next frame. */
        int temporal_diff2 = ( ABS(next[-refs] - cur[-refs]) + ABS(next[+refs] - cur[+refs]) ) >> 1;
        /* For the actual difference, use the largest of the previous average diffs. */
        int diff           = MAX3(temporal_diff0>>1, temporal_diff1, temporal_diff2);

        int spatial_pred;
        
        if( eedi2_mode )
        {
            /* Who needs yadif's spatial predictions when we can have EEDI2's? */
            spatial_pred = eedi2_guess[0];
            eedi2_guess++;
        }
        else // Yadif spatial interpolation
        {
            /* SAD of how the pixel-1, the pixel, and the pixel+1 change from the line above to below. */ 
            int spatial_score  = ABS(cur[-refs-1] - cur[+refs-1]) + ABS(cur[-refs]-cur[+refs]) +
                                         ABS(cur[-refs+1] - cur[+refs+1]) - 1;         
            
            /* Spatial pred is either a bilinear or cubic vertical interpolation. */
            if( ( pv->mode & MODE_CUBIC ) && !vertical_edge)
            {
                spatial_pred = cubic_interpolate_pixel( cur[-3*refs], cur[-refs], cur[+refs], cur[3*refs] );
            }
            else
            {
                spatial_pred = (c+e)>>1;
            }

        /* EDDI: Edge Directed Deinterlacing Interpolation
           Checks 4 different slopes to see if there is more similarity along a diagonal
           than there was vertically. If a diagonal is more similar, then it indicates
           an edge, so interpolate along that instead of a vertical line, using either
           linear or cubic interpolation depending on mode. */
        #define YADIF_CHECK(j)\
                {   int score = ABS(cur[-refs-1+j] - cur[+refs-1-j])\
                              + ABS(cur[-refs  +j] - cur[+refs  -j])\
                              + ABS(cur[-refs+1+j] - cur[+refs+1-j]);\
                    if( score < spatial_score ){\
                        spatial_score = score;\
                        if( ( pv->mode & MODE_CUBIC ) && !vertical_edge )\
                        {\
                            switch(j)\
                            {\
                                case -1:\
                                    spatial_pred = cubic_interpolate_pixel(cur[-3 * refs - 3], cur[-refs -1], cur[+refs + 1], cur[3* refs + 3] );\
                                break;\
                                case -2:\
                                    spatial_pred = cubic_interpolate_pixel( ( ( cur[-3*refs - 4] + cur[-refs - 4] ) / 2 ) , cur[-refs -2], cur[+refs + 2], ( ( cur[3*refs + 4] + cur[refs + 4] ) / 2 ) );\
                                break;\
                                case 1:\
                                    spatial_pred = cubic_interpolate_pixel(cur[-3 * refs +3], cur[-refs +1], cur[+refs - 1], cur[3* refs -3] );\
                                break;\
                                case 2:\
                                    spatial_pred = cubic_interpolate_pixel(( ( cur[-3*refs + 4] + cur[-refs + 4] ) / 2 ), cur[-refs +2], cur[+refs - 2], ( ( cur[3*refs - 4] + cur[refs - 4] ) / 2 ) );\
                                break;\
                            }\
                        }\
                        else\
                        {\
                            spatial_pred = ( cur[-refs +j] + cur[+refs -j] ) >>1;\
                        }\

                        if( x >= 2 && x <= w - 3 )
                        {
                            YADIF_CHECK(-1)
                            if( x >= 3 && x <= w - 4 )
                            {
                                YADIF_CHECK(-2) }} }}
                            }
                        }
                        if( x >= 2 && x <= w - 3 )
                        {
                            YADIF_CHECK(1)
                            if( x >= 3 && x <= w - 4 )
                            {
                                YADIF_CHECK(2) }} }}
                            }
                        }
        }

        /* Temporally adjust the spatial prediction by
           comparing against lines in the adjacent fields. */
        int b = (prev2[-2*refs] + next2[-2*refs])>>1;
        int f = (prev2[+2*refs] + next2[+2*refs])>>1;
        
        /* Find the median value */
        int max = MAX3(d-e, d-c, MIN(b-c, f-e));
        int min = MIN3(d-e, d-c, MAX(b-c, f-e));
        diff = MAX3( diff, min, -max );
        
        if( spatial_pred > d + diff )
        {
            spatial_pred = d + diff;
        }
        else if( spatial_pred < d - diff )
        {
            spatial_pred = d - diff;
        }
        
        dst[0] = spatial_pred;
                        
        dst++;
        cur++;
        prev++;
        next++;
        prev2++;
        next2++;
    }
}

/*
 * deinterlace this segment of all three planes in a single thread.
 */
void yadif_decomb_filter_thread( void *thread_args_v )
{
    yadif_arguments_t *yadif_work = NULL;
    hb_filter_private_t * pv;
    int run = 1;
    int plane;
    int segment, segment_start, segment_stop;
    yadif_thread_arg_t *thread_args = thread_args_v;
    uint8_t **dst;
    int parity, tff, y, w, h, penultimate, ultimate, ref_stride, is_combed;

    pv = thread_args->pv;
    segment = thread_args->segment;

    hb_log("yadif thread started for segment %d", segment);

    while( run )
    {
        /*
         * Wait here until there is work to do. hb_lock() blocks until
         * render releases it to say that there is more work to do.
         */
        hb_lock( pv->yadif_begin_lock[segment] );

        yadif_work = &pv->yadif_arguments[segment];

        if( yadif_work->stop )
        {
            /*
             * No more work to do, exit this thread.
             */
            run = 0;
            continue;
        } 

        if( yadif_work->dst == NULL )
        {
            hb_error( "thread started when no work available" );
            hb_snooze(500);
            continue;
        }
        
        is_combed = pv->yadif_arguments[segment].is_combed;

        /*
         * Process all three planes, but only this segment of it.
         */
        for( plane = 0; plane < 3; plane++)
        {

            dst = yadif_work->dst;
            parity = yadif_work->parity;
            tff = yadif_work->tff;
            w = pv->width[plane];
            h = pv->height[plane];
            penultimate = h - 2;
            ultimate = h - 1;
            ref_stride = pv->ref_stride[plane];
            segment_start = ( h / pv->cpu_count ) * segment;
            int size_of_w = w * sizeof(uint8_t);
            
            if( segment == pv->cpu_count - 1 )
            {
                /*
                 * Final segment
                 */
                segment_stop = h;
            } else {
                segment_stop = ( h / pv->cpu_count ) * ( segment + 1 );
            }

            for( y = segment_start; y < segment_stop; y++ )
            {
                if( is_combed == 2 )
                {
                    /* This line gets blend filtered, not yadif filtered. */
                    uint8_t *cur  = &pv->ref[1][plane][y*ref_stride];
                    uint8_t *dst2 = &dst[plane][y*w];
                    /* These will be useful if we ever do temporal blending. */
                    // uint8_t *prev = &pv->ref[0][plane][y*ref_stride];
                    // uint8_t *next = &pv->ref[2][plane][y*ref_stride];

                    blend_filter_line( dst2, cur, plane, y, pv );
                }
                else if( pv->mode == MODE_CUBIC && is_combed && ( ( y ^ parity ) & 1 ) )
                {
                    /* Just apply vertical cubic interpolation */
                    uint8_t *cur  = &pv->ref[1][plane][y*ref_stride];
                    uint8_t *dst2 = &dst[plane][y*w];
                    
                    cubic_interpolate_line( dst2, cur, plane, y, pv );
                }
                else if( pv->mode & MODE_YADIF && ( ( y ^ parity ) &  1 )  && ( is_combed == 1 ) )
                {
                    /* This line gets yadif filtered. It is the bottom field
                       when TFF and vice-versa. It's the field that gets
                       filtered. Because yadif needs 2 lines above and below
                       the one being filtered, we need to mirror the edges.
                       When TFF, this means replacing the 2nd line with a
                       copy of the 1st, and the last with the second-to-last. */
                    if( y > 1 && y < ( h -2 ) )
                    {
                        /* This isn't the top or bottom, proceed as normal to yadif. */
                        uint8_t *prev = &pv->ref[0][plane][y*ref_stride];
                        uint8_t *cur  = &pv->ref[1][plane][y*ref_stride];
                        uint8_t *next = &pv->ref[2][plane][y*ref_stride];
                        uint8_t *dst2 = &dst[plane][y*w];

                        yadif_filter_line( dst2, 
                                           prev, 
                                           cur, 
                                           next, 
                                           plane, 
                                           parity ^ tff,
                                           y, 
                                           pv );
                    }
                    else if( y == 0 )
                    {
                        /* BFF, so y0 = y1 */
                        memcpy( &dst[plane][y*w],
                                &pv->ref[1][plane][1*ref_stride],
                                size_of_w );
                    }
                    else if( y == 1 )
                    {
                        /* TFF, so y1 = y0 */
                        memcpy( &dst[plane][y*w],
                                &pv->ref[1][plane][0],
                                size_of_w );
                    }
                    else if( y == penultimate )
                    {
                        /* BFF, so penultimate y = ultimate y */
                        memcpy( &dst[plane][y*w],
                                &pv->ref[1][plane][ultimate*ref_stride],
                                size_of_w );
                    }
                    else if( y == ultimate )
                    {
                        /* TFF, so ultimate y = penultimate y */
                        memcpy( &dst[plane][y*w],
                                &pv->ref[1][plane][penultimate*ref_stride],
                                size_of_w );
                    }
                }
                else
                {
                    memcpy( &dst[plane][y*w],
                            &pv->ref[1][plane][y*ref_stride],
                            size_of_w );              
                }
            }
        }
        /*
         * Finished this segment, let everyone know.
         */
        hb_unlock( pv->yadif_complete_lock[segment] );
    }
    free( thread_args_v );
}

static void yadif_filter( uint8_t ** dst,
                          int parity,
                          int tff,
                          hb_filter_private_t * pv )
{
    /* If we're running comb detection, do it now, otherwise default to true. */
    int is_combed = pv->spatial_metric >= 0 ? comb_segmenter( pv ) : 1;
    
    /* The comb detector suggests three different values:
       0: Don't comb this frame.
       1: Deinterlace this frame.
       2: Blend this frame.
       Since that might conflict with the filter's mode,
       it may be necesary to adjust this value.          */
    if( is_combed == 1 && (pv->mode == MODE_BLEND) )
    {
        /* All combed frames are getting blended */
        is_combed = 2;
    }
    else if( is_combed == 2 && !( pv->mode & MODE_BLEND ) )
    {
        /* Blending is disabled, so force interpolation of these frames. */
        is_combed = 1;
    }
    if( is_combed == 1 &&
        ( pv->mode & MODE_BLEND ) &&
        !( pv->mode & ( MODE_YADIF | MODE_EEDI2 | MODE_CUBIC ) ) )
    {
        /* Deinterlacers are disabled, blending isn't, so blend these frames. */
        is_combed = 2;
    }
    else if( is_combed &&
             !( pv->mode & ( MODE_BLEND | MODE_YADIF | MODE_EEDI2 | MODE_CUBIC | MODE_MASK ) ) )
    {
        /* No deinterlacer or mask chosen, pass the frame through. */
        is_combed = 0;
    }
    
    if( is_combed == 1 )
    {
        pv->deinterlaced_frames++;
    }
    else if( is_combed == 2 )
    {
        pv->blended_frames++;
    }
    else
    {
        pv->unfiltered_frames++;
    }
    
    if( is_combed == 1 && ( pv->mode & MODE_EEDI2 ) )
    {
        /* Generate an EEDI2 interpolation */
        eedi2_planer( pv );
    }
    
    if( is_combed )
    {
        if( ( pv->mode & MODE_EEDI2 ) && !( pv->mode & MODE_YADIF ) && is_combed == 1 )
        {
            // Just pass through the EEDI2 interpolation
            int i;
            for( i = 0; i < 3; i++ )
            {
                uint8_t * ref = pv->eedi_full[DST2PF][i];
                uint8_t * dest = dst[i];

                int w = pv->width[i];
                int ref_stride = pv->ref_stride[i];

                int y;
                for( y = 0; y < pv->height[i]; y++ )
                {
                    memcpy(dest, ref, w);
                    dest += w;
                    ref += ref_stride;
                }
            }
        }
        else
        {
            int segment;

            for( segment = 0; segment < pv->cpu_count; segment++ )
            {  
                /*
                 * Setup the work for this plane.
                 */
                pv->yadif_arguments[segment].parity = parity;
                pv->yadif_arguments[segment].tff = tff;
                pv->yadif_arguments[segment].dst = dst;
                pv->yadif_arguments[segment].is_combed = is_combed;

                /*
                 * Let the thread for this plane know that we've setup work 
                 * for it by releasing the begin lock (ensuring that the
                 * complete lock is already locked so that we block when
                 * we try to lock it again below).
                 */
                hb_lock( pv->yadif_complete_lock[segment] );
                hb_unlock( pv->yadif_begin_lock[segment] );
            }

            /*
             * Wait until all three threads have completed by trying to get
             * the complete lock that we locked earlier for each thread, which
             * will block until that thread has completed the work on that
             * plane.
             */
            for( segment = 0; segment < pv->cpu_count; segment++ )
            {
                hb_lock( pv->yadif_complete_lock[segment] );
                hb_unlock( pv->yadif_complete_lock[segment] );
            }

            /*
             * Entire frame is now deinterlaced.
             */
        }
    }
    else
    {
        /*  Just passing through... */
        
        /* For mcdeint's benefit... */
        pv->yadif_arguments[0].is_combed = is_combed; // 0
        
        int i;
        for( i = 0; i < 3; i++ )
        {
            uint8_t * ref = pv->ref[1][i];
            uint8_t * dest = dst[i];
            
            int w = pv->width[i];
            int h = pv->height[i];
            int ref_stride = pv->ref_stride[i];
            
            int y;
            for( y = 0; y < h; y++ )
            {
                memcpy(dest, ref, w);
                dest += w;
                ref += ref_stride;
            }
        }
    }
    
    if( pv->mode & MODE_MASK && pv->spatial_metric >= 0 )
    {
        if( pv->mode == MODE_MASK || (pv->mode & MODE_MASK && pv->mode & MODE_GAMMA) || is_combed || (pv->mode & MODE_MASK && pv->mode & MODE_FILTER ))
        apply_mask( pv );
    }
}

static void mcdeint_filter( uint8_t ** dst,
                            uint8_t ** src,
                            int parity,
                            hb_filter_private_t * pv )
{
    int x, y, i;
    int out_size;

#ifdef SUPPRESS_AV_LOG
    /* TODO: temporarily change log level to suppress obnoxious debug output */
    int loglevel = av_log_get_level();
    av_log_set_level( AV_LOG_QUIET );
#endif

    for( i=0; i<3; i++ )
    {
        pv->mcdeint_frame->data[i] = src[i];
        pv->mcdeint_frame->linesize[i] = pv->width[i];
    }
    pv->mcdeint_avctx_enc->me_cmp     = FF_CMP_SAD;
    pv->mcdeint_avctx_enc->me_sub_cmp = FF_CMP_SAD;
    pv->mcdeint_frame->quality        = pv->mcdeint_qp * FF_QP2LAMBDA;

    out_size = avcodec_encode_video( pv->mcdeint_avctx_enc,
                                     pv->mcdeint_outbuf,
                                     pv->mcdeint_outbuf_size,
                                     pv->mcdeint_frame );

    pv->mcdeint_frame_dec = pv->mcdeint_avctx_enc->coded_frame;

    for( i = 0; i < 3; i++ )
    {
        int w    = pv->width[i];
        int h    = pv->height[i];
        int fils = pv->mcdeint_frame_dec->linesize[i];
        int srcs = pv->width[i];

        for( y = 0; y < h; y++ )
        {
            if( (y ^ parity) & 1 )
            {
                for( x = 0; x < w; x++ )
                {
                    if( (x-1)+(y-1)*w >= 0 && (x+1)+(y+1)*w < w*h )
                    {
                        uint8_t * filp =
                            &pv->mcdeint_frame_dec->data[i][x + y*fils];
                        uint8_t * srcp = &src[i][x + y*srcs];

                        int diff0 = filp[-fils] - srcp[-srcs];
                        int diff1 = filp[+fils] - srcp[+srcs];
                        int spatial_score;
                        
                        spatial_score =
                            ABS(srcp[-srcs-1] - srcp[+srcs-1]) +
                            ABS(srcp[-srcs  ] - srcp[+srcs  ]) +
                            ABS(srcp[-srcs+1] - srcp[+srcs+1]) - 1;

                        int temp = filp[0];

#define MCDEINT_CHECK(j)\
                        {   int score = ABS(srcp[-srcs-1+j] - srcp[+srcs-1-j])\
                                      + ABS(srcp[-srcs  +j] - srcp[+srcs  -j])\
                                      + ABS(srcp[-srcs+1+j] - srcp[+srcs+1-j]);\
                            if( score < spatial_score ) {\
                                spatial_score = score;\
                                diff0 = filp[-fils+j] - srcp[-srcs+j];\
                                diff1 = filp[+fils-j] - srcp[+srcs-j];

                        if( x >= 2 && x <= w - 3 )
                        {
                            MCDEINT_CHECK(-1)
                            if( x >= 3 && x <= w - 4 )
                            {
                                MCDEINT_CHECK(-2) }} }}
                            }
                        }
                        if( x >= 2 && x <= w - 3 )
                        {
                            MCDEINT_CHECK(1)
                            if( x >= 3 && x <= w - 4 )
                            {
                                MCDEINT_CHECK(2) }} }}
                            }
                        }

                        if(diff0 + diff1 > 0)
                        {
                            temp -= (diff0 + diff1 -
                                     ABS( ABS(diff0) - ABS(diff1) ) / 2) / 2;
                        }
                        else
                        {
                            temp -= (diff0 + diff1 +
                                     ABS( ABS(diff0) - ABS(diff1) ) / 2) / 2;
                        }

                        filp[0] = dst[i][x + y*w] =
                            temp > 255U ? ~(temp>>31) : temp;
                    }
                    else
                    {
                        dst[i][x + y*w] =
                            pv->mcdeint_frame_dec->data[i][x + y*fils];
                    }
                }
            }
            else
            {
                for( x = 0; x < w; x++ )
                {
                    pv->mcdeint_frame_dec->data[i][x + y*fils] =
                        dst[i][x + y*w]= src[i][x + y*srcs];
                }
            }
        }
    }

#ifdef SUPPRESS_AV_LOG
    /* TODO: restore previous log level */
    av_log_set_level(loglevel);
#endif
}

hb_filter_private_t * hb_decomb_init( int pix_fmt,
                                           int width,
                                           int height,
                                           char * settings )
{
    if( pix_fmt != PIX_FMT_YUV420P )
    {
        return 0;
    }

    hb_filter_private_t * pv = calloc( 1, sizeof(struct hb_filter_private_s) );

    pv->pix_fmt = pix_fmt;

    pv->width[0]  = width;
    pv->height[0] = height;
    pv->width[1]  = pv->width[2]  = width >> 1;
    pv->height[1] = pv->height[2] = height >> 1;

    build_gamma_lut( pv );
    
    pv->buf_out[0] = hb_video_buffer_init( width, height );
    pv->buf_out[1] = hb_video_buffer_init( width, height );
    pv->buf_settings = hb_buffer_init( 0 );

    pv->deinterlaced_frames = 0;
    pv->blended_frames = 0;
    pv->unfiltered_frames = 0;

    pv->yadif_ready    = 0;

    pv->mode     = MODE_YADIF | MODE_BLEND | MODE_CUBIC;
    pv->filter_mode = FILTER_ERODE_DILATE;
    pv->spatial_metric = 2;
    pv->motion_threshold = 6;
    pv->spatial_threshold = 9;
    pv->block_threshold = 80;
    pv->block_width = 16;
    pv->block_height = 16;
    
    pv->magnitude_threshold = 10;
    pv->variance_threshold = 20;
    pv->laplacian_threshold = 20;
    pv->dilation_threshold = 4;
    pv->erosion_threshold = 2;
    pv->noise_threshold = 50;
    pv->maximum_search_distance = 24;
    pv->post_processing = 1;

    pv->parity   = PARITY_DEFAULT;

    pv->mcdeint_mode   = MCDEINT_MODE_DEFAULT;
    pv->mcdeint_qp     = MCDEINT_QP_DEFAULT;

    if( settings )
    {
        sscanf( settings, "%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d:%d",
                &pv->mode,
                &pv->spatial_metric,
                &pv->motion_threshold,
                &pv->spatial_threshold,
                &pv->filter_mode,
                &pv->block_threshold,
                &pv->block_width,
                &pv->block_height,
                &pv->magnitude_threshold,
                &pv->variance_threshold,
                &pv->laplacian_threshold,
                &pv->dilation_threshold,
                &pv->erosion_threshold,
                &pv->noise_threshold,
                &pv->maximum_search_distance,
                &pv->post_processing,
                &pv->parity );
    }
    
    pv->cpu_count = hb_get_cpu_count();
    

    if( pv->mode & MODE_MCDEINT )
    {
        pv->mcdeint_mode = 2;
    }
    
    /* Allocate yadif specific buffers */
    int i, j;
    for( i = 0; i < 3; i++ )
    {
        int is_chroma = !!i;
        int w = ((width   + 31) & (~31))>>is_chroma;
        int h = ((height+6+ 31) & (~31))>>is_chroma;

        pv->ref_stride[i] = w;

        for( j = 0; j < 3; j++ )
        {
            pv->ref[j][i] = calloc( 1, w*h*sizeof(uint8_t) ) + 3*w;
        }
    }

    /* Allocate buffers to store comb masks. */
    for( i = 0; i < 3; i++ )
    {
        int is_chroma = !!i;
        int w = ((pv->width[0]   + 31) & (~31))>>is_chroma;
        int h = ((pv->height[0]+6+ 31) & (~31))>>is_chroma;

        pv->mask[i] = calloc( 1, w*h*sizeof(uint8_t) ) + 3*w;
        pv->mask_filtered[i] = calloc( 1, w*h*sizeof(uint8_t) ) + 3*w;
        pv->mask_temp[i] = calloc( 1, w*h*sizeof(uint8_t) ) + 3*w;
    }
    
    if( pv->mode & MODE_EEDI2 )
    {
        /* Allocate half-height eedi2 buffers */
        height = pv->height[0] / 2;
        for( i = 0; i < 3; i++ )
        {
            int is_chroma = !!i;
            int w = ((width   + 31) & (~31))>>is_chroma;
            int h = ((height+6+ 31) & (~31))>>is_chroma;

            for( j = 0; j < 4; j++ )
            {
                pv->eedi_half[j][i] = calloc( 1, w*h*sizeof(uint8_t) ) + 3*w;
            }
        }

        /* Allocate full-height eedi2 buffers */
        height = pv->height[0];
        for( i = 0; i < 3; i++ )
        {
            int is_chroma = !!i;
            int w = ((width   + 31) & (~31))>>is_chroma;
            int h = ((height+6+ 31) & (~31))>>is_chroma;

            for( j = 0; j < 5; j++ )
            {
                pv->eedi_full[j][i] = calloc( 1, w*h*sizeof(uint8_t) ) + 3*w;
            }
        }
    }
    
     /*
      * Create yadif threads and locks.
      */
     pv->yadif_threads = malloc( sizeof( hb_thread_t* ) * pv->cpu_count );
     pv->yadif_begin_lock = malloc( sizeof( hb_lock_t * ) * pv->cpu_count );
     pv->yadif_complete_lock = malloc( sizeof( hb_lock_t * ) * pv->cpu_count );
     pv->yadif_arguments = malloc( sizeof( yadif_arguments_t ) * pv->cpu_count );

     for( i = 0; i < pv->cpu_count; i++ )
     {
         yadif_thread_arg_t *thread_args;

         thread_args = malloc( sizeof( yadif_thread_arg_t ) );

         if( thread_args )
         {
             thread_args->pv = pv;
             thread_args->segment = i;

             pv->yadif_begin_lock[i] = hb_lock_init();
             pv->yadif_complete_lock[i] = hb_lock_init();

             /*
              * Important to start off with the threads locked waiting
              * on input.
              */
             hb_lock( pv->yadif_begin_lock[i] );

             pv->yadif_arguments[i].stop = 0;
             pv->yadif_arguments[i].dst = NULL;
             
             pv->yadif_threads[i] = hb_thread_init( "yadif_filter_segment",
                                                    yadif_decomb_filter_thread,
                                                    thread_args,
                                                    HB_NORMAL_PRIORITY );
         }
         else
         {
             hb_error( "yadif could not create threads" );
         }
    }
    
    /*
     * Create decomb threads and locks.
     */
    pv->decomb_threads = malloc( sizeof( hb_thread_t* ) * pv->cpu_count );
    pv->decomb_begin_lock = malloc( sizeof( hb_lock_t * ) * pv->cpu_count );
    pv->decomb_complete_lock = malloc( sizeof( hb_lock_t * ) * pv->cpu_count );
    pv->decomb_arguments = malloc( sizeof( decomb_arguments_t ) * pv->cpu_count );
    
    for( i = 0; i < pv->cpu_count; i++ )
    {
        decomb_thread_arg_t *decomb_thread_args;
    
        decomb_thread_args = malloc( sizeof( decomb_thread_arg_t ) );
    
        if( decomb_thread_args )
        {
            decomb_thread_args->pv = pv;
            decomb_thread_args->segment = i;
    
            pv->decomb_begin_lock[i] = hb_lock_init();
            pv->decomb_complete_lock[i] = hb_lock_init();
    
            /*
             * Important to start off with the threads locked waiting
             * on input.
             */
            hb_lock( pv->decomb_begin_lock[i] );
    
            pv->decomb_arguments[i].stop = 0;
    
            pv->decomb_threads[i] = hb_thread_init( "decomb_filter_segment",
                                                   decomb_filter_thread,
                                                   decomb_thread_args,
                                                   HB_NORMAL_PRIORITY );
        }
        else
        {
            hb_error( "decomb could not create threads" );
        }
    }
    
    if( pv->mode & MODE_EEDI2 )
    {
        /*
         * Create eedi2 threads and locks.
         */
        pv->eedi2_threads = malloc( sizeof( hb_thread_t* ) * 3 );
        pv->eedi2_begin_lock = malloc( sizeof( hb_lock_t * ) * 3 );
        pv->eedi2_complete_lock = malloc( sizeof( hb_lock_t * ) * 3 );
        pv->eedi2_arguments = malloc( sizeof( eedi2_arguments_t ) * 3 );

        if( pv->post_processing > 1 )
        {
            pv->cx2 = (int*)eedi2_aligned_malloc(pv->height[0]*pv->ref_stride[0]*sizeof(int), 16);
            pv->cy2 = (int*)eedi2_aligned_malloc(pv->height[0]*pv->ref_stride[0]*sizeof(int), 16);
            pv->cxy = (int*)eedi2_aligned_malloc(pv->height[0]*pv->ref_stride[0]*sizeof(int), 16);
            pv->tmpc = (int*)eedi2_aligned_malloc(pv->height[0]*pv->ref_stride[0]*sizeof(int), 16);
            if( !pv->cx2 || !pv->cy2 || !pv->cxy || !pv->tmpc )
                hb_log("EEDI2: failed to malloc derivative arrays");
            else
                hb_log("EEDI2: successfully mallloced derivative arrays");
        }

        for( i = 0; i < 3; i++ )
        {
            eedi2_thread_arg_t *eedi2_thread_args;

            eedi2_thread_args = malloc( sizeof( eedi2_thread_arg_t ) );

            if( eedi2_thread_args )
            {
                eedi2_thread_args->pv = pv;
                eedi2_thread_args->plane = i;

                pv->eedi2_begin_lock[i] = hb_lock_init();
                pv->eedi2_complete_lock[i] = hb_lock_init();

                /*
                 * Important to start off with the threads locked waiting
                 * on input.
                 */
                hb_lock( pv->eedi2_begin_lock[i] );

                pv->eedi2_arguments[i].stop = 0;

                pv->eedi2_threads[i] = hb_thread_init( "eedi2_filter_segment",
                                                       eedi2_filter_thread,
                                                       eedi2_thread_args,
                                                       HB_NORMAL_PRIORITY );
            }
            else
            {
                hb_error( "eedi2 could not create threads" );
            }
        }
    }
    
    
    /* Allocate mcdeint specific buffers */
    if( pv->mcdeint_mode >= 0 )
    {
        avcodec_init();
        avcodec_register_all();
        AVCodec * enc = avcodec_find_encoder( CODEC_ID_SNOW );
        int i;
        for (i = 0; i < 3; i++ )
        {
            AVCodecContext * avctx_enc;

            avctx_enc = pv->mcdeint_avctx_enc = avcodec_alloc_context();

            avctx_enc->width                    = width;
            avctx_enc->height                   = height;
            avctx_enc->time_base                = (AVRational){1,25};  // meaningless
            avctx_enc->gop_size                 = 300;
            avctx_enc->max_b_frames             = 0;
            avctx_enc->pix_fmt                  = PIX_FMT_YUV420P;
            avctx_enc->flags                    = CODEC_FLAG_QSCALE | CODEC_FLAG_LOW_DELAY;
            avctx_enc->strict_std_compliance    = FF_COMPLIANCE_EXPERIMENTAL;
            avctx_enc->global_quality           = 1;
            avctx_enc->flags2                   = CODEC_FLAG2_MEMC_ONLY;
            avctx_enc->me_cmp                   = FF_CMP_SAD; //SSE;
            avctx_enc->me_sub_cmp               = FF_CMP_SAD; //SSE;
            avctx_enc->mb_cmp                   = FF_CMP_SSE;

            switch( pv->mcdeint_mode )
            {
                case 3:
                    avctx_enc->refs = 3;
                case 2:
                    avctx_enc->me_method = ME_ITER;
                case 1:
                    avctx_enc->flags |= CODEC_FLAG_4MV;
                    avctx_enc->dia_size =2;
                case 0:
                    avctx_enc->flags |= CODEC_FLAG_QPEL;
            }

            hb_avcodec_open(avctx_enc, enc, 0);
        }

        pv->mcdeint_frame       = avcodec_alloc_frame();
        pv->mcdeint_outbuf_size = width * height * 10;
        pv->mcdeint_outbuf      = malloc( pv->mcdeint_outbuf_size );
    }

    return pv;
}

void hb_decomb_close( hb_filter_private_t * pv )
{
    if( !pv )
    {
        return;
    }
    
    hb_log("decomb: deinterlaced %i | blended %i | unfiltered %i | total %i", pv->deinterlaced_frames, pv->blended_frames, pv->unfiltered_frames, pv->deinterlaced_frames + pv->blended_frames + pv->unfiltered_frames);

    /* Cleanup frame buffers */
    if( pv->buf_out[0] )
    {
        hb_buffer_close( &pv->buf_out[0] );
    }
    if( pv->buf_out[1] )
    {
        hb_buffer_close( &pv->buf_out[1] );
    }
    if (pv->buf_settings )
    {
        hb_buffer_close( &pv->buf_settings );
    }

    /* Cleanup yadif specific buffers */
    int i;
    for( i = 0; i<3*3; i++ )
    {
        uint8_t **p = &pv->ref[i%3][i/3];
        if (*p)
        {
            free( *p - 3*pv->ref_stride[i/3] );
            *p = NULL;
        }
    }
    
    /* Cleanup combing masks. */
    for( i = 0; i<3*3; i++ )
    {
        uint8_t **p = &pv->mask[i/3];
        uint8_t **p2 = &pv->mask_filtered[i/3];
        uint8_t **p3 = &pv->mask_temp[i/3];

        if (*p)
        {
            free( *p - 3*pv->ref_stride[i/3] );
            *p = NULL;
        }
        if (*p2)
        {
            free( *p2 - 3*pv->ref_stride[i/3] );
            *p2 = NULL;
        }
        if (*p3)
        {
            free( *p3 - 3*pv->ref_stride[i/3] );
            *p3 = NULL;
        }
    }
    
    if( pv->mode & MODE_EEDI2 )
    {
        /* Cleanup eedi-half  buffers */
        int j;
        for( i = 0; i<3; i++ )
        {
            for( j = 0; j < 4; j++ )
            {
                uint8_t **p = &pv->eedi_half[j][i];
                if (*p)
                {
                    free( *p - 3*pv->ref_stride[i] );
                    *p = NULL;
                }            
            }
        }

        /* Cleanup eedi-full  buffers */
        for( i = 0; i<3; i++ )
        {
            for( j = 0; j < 5; j++ )
            {
                uint8_t **p = &pv->eedi_full[j][i];
                if (*p)
                {
                    free( *p - 3*pv->ref_stride[i] );
                    *p = NULL;
                }            
            }
        }
    }
    
    if( pv->post_processing > 1  && ( pv->mode & MODE_EEDI2 ) )
    {
        if (pv->cx2) eedi2_aligned_free(pv->cx2);
        if (pv->cy2) eedi2_aligned_free(pv->cy2);
        if (pv->cxy) eedi2_aligned_free(pv->cxy);
        if (pv->tmpc) eedi2_aligned_free(pv->tmpc);
    }
    
    for( i = 0; i < pv->cpu_count; i++)
    {
        /*
         * Tell each yadif thread to stop, and then cleanup.
         */
        pv->yadif_arguments[i].stop = 1;
        hb_unlock(  pv->yadif_begin_lock[i] );

        hb_thread_close( &pv->yadif_threads[i] );
        hb_lock_close( &pv->yadif_begin_lock[i] );
        hb_lock_close( &pv->yadif_complete_lock[i] );
    }
    
    /*
     * free memory for yadif structs
     */
    free( pv->yadif_threads );
    free( pv->yadif_begin_lock );
    free( pv->yadif_complete_lock );
    free( pv->yadif_arguments );
    
    for( i = 0; i < pv->cpu_count; i++)
    {
        /*
         * Tell each decomb thread to stop, and then cleanup.
         */
        pv->decomb_arguments[i].stop = 1;
        hb_unlock(  pv->decomb_begin_lock[i] );

        hb_thread_close( &pv->decomb_threads[i] );
        hb_lock_close( &pv->decomb_begin_lock[i] );
        hb_lock_close( &pv->decomb_complete_lock[i] );
    }
    
    /*
     * free memory for decomb structs
     */
    free( pv->decomb_threads );
    free( pv->decomb_begin_lock );
    free( pv->decomb_complete_lock );
    free( pv->decomb_arguments );
    
    if( pv->mode & MODE_EEDI2 )
    {
        for( i = 0; i < 3; i++)
        {
            /*
             * Tell each eedi2 thread to stop, and then cleanup.
             */
            pv->eedi2_arguments[i].stop = 1;
            hb_unlock(  pv->eedi2_begin_lock[i] );

            hb_thread_close( &pv->eedi2_threads[i] );
            hb_lock_close( &pv->eedi2_begin_lock[i] );
            hb_lock_close( &pv->eedi2_complete_lock[i] );
        }

        /*
         * free memory for eedi2 structs
         */
        free( pv->eedi2_threads );
        free( pv->eedi2_begin_lock );
        free( pv->eedi2_complete_lock );
        free( pv->eedi2_arguments );
    }
    
    /* Cleanup mcdeint specific buffers */
    if( pv->mcdeint_mode >= 0 )
    {
        if( pv->mcdeint_avctx_enc )
        {
            hb_avcodec_close( pv->mcdeint_avctx_enc );
            av_freep( &pv->mcdeint_avctx_enc );
        }
        if( pv->mcdeint_outbuf )
        {
            free( pv->mcdeint_outbuf );
        }
    }

    free( pv );
}

int hb_decomb_work( const hb_buffer_t * cbuf_in,
                    hb_buffer_t ** buf_out,
                    int pix_fmt,
                    int width,
                    int height,
                    hb_filter_private_t * pv )
{
    hb_buffer_t * buf_in = (hb_buffer_t *)cbuf_in;

    if( !pv ||
        pix_fmt != pv->pix_fmt ||
        width   != pv->width[0] ||
        height  != pv->height[0] )
    {
        return FILTER_FAILED;
    }

    avpicture_fill( &pv->pic_in, buf_in->data,
                    pix_fmt, width, height );

    /* Determine if top-field first layout */
    int tff;
    if( pv->parity < 0 )
    {
        tff = !!(buf_in->flags & PIC_FLAG_TOP_FIELD_FIRST);
    }
    else
    {
        tff = (pv->parity & 1) ^ 1;
    }

    /* Store current frame in yadif cache */
    store_ref( (const uint8_t**)pv->pic_in.data, pv );

    /* If yadif is not ready, store another ref and return FILTER_DELAY */
    if( pv->yadif_ready == 0 )
    {
        store_ref( (const uint8_t**)pv->pic_in.data, pv );

        hb_buffer_copy_settings( pv->buf_settings, buf_in );

        /* don't let 'work_loop' send a chapter mark upstream */
        buf_in->new_chap  = 0;

        pv->yadif_ready = 1;

        return FILTER_DELAY;
    }

    /* Perform yadif filtering */        
    int frame;
    for( frame = 0; frame <= ( ( pv->mode & MODE_MCDEINT ) ? 1 : 0 ) ; frame++ )
// This would be what to use for bobbing: for( frame = 0; frame <= 0 ; frame++ )
    {

#if 0        
        /* Perhaps skip the second run if the frame is uncombed? */
        if( frame && !pv->yadif_arguments[0].is_combed )
        {
            break;
        }
#endif        
        int parity = frame ^ tff ^ 1;

// This will be for bobbing
#if 0
        if( pv->alternator )
        {
            parity = !parity;
            pv->alternator = 0;
        }
        else
        {
            pv->alternator = 1;
        }
#endif
        pv->tff = !parity;

        avpicture_fill( &pv->pic_out, pv->buf_out[!(frame^1)]->data,
                        pix_fmt, width, height );

        /* XXX
            Should check here and only bother filtering field 2 when
           field 1 was detected as combed.
           And when it's not, it's a progressive frame,
           so mcdeint should be skipped...
        */
        yadif_filter( pv->pic_out.data, parity, tff, pv );

        /* Commented out code in the line below would skip mcdeint
           on uncombed frames. Possibly a bad idea, since mcdeint
           maintains the same snow context for the entire video... */
        if( pv->mcdeint_mode >= 0 /* && pv->yadif_arguments[0].is_combed */)
        {
            /* Perform mcdeint filtering */
            avpicture_fill( &pv->pic_in,  pv->buf_out[(frame^1)]->data,
                            pix_fmt, width, height );

            mcdeint_filter( pv->pic_in.data, pv->pic_out.data, parity, pv );
        }

        *buf_out = pv->buf_out[!(frame^1)];
    }

    /* Copy buffered settings to output buffer settings */
    hb_buffer_copy_settings( *buf_out, pv->buf_settings );

    /* Replace buffered settings with input buffer settings */
    hb_buffer_copy_settings( pv->buf_settings, buf_in );

    /* don't let 'work_loop' send a chapter mark upstream */
    buf_in->new_chap  = 0;

    return FILTER_OK;
}
