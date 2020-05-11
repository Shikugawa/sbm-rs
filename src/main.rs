use rand::Rng;
use std::num;

struct CompleteGraph {
    node_num: u32,
    weights: Vec<Vec<f32>>,
}

impl CompleteGraph {
    fn new(node_num: u32) -> Self {
        let mut rng = rand::thread_rng();
        let mut weights: Vec<Vec<f32>> = Vec::with_capacity(node_num as usize);

        for i in 0..node_num {
            let mut col: Vec<f32> = Vec::with_capacity(node_num as usize);
            for j in 0..node_num {
                col.push(0.0);
            }
            weights.push(col);
        }

        for i in 0..node_num {
            for j in 0..node_num {
                if weights[i as usize][j as usize] != 0
                    || weights[i as usize][j as usize] == weights[j as usize][i as usize]
                {
                    continue;
                }
                let weight: f32 = rng.gen_range(1.0, 10000.0);
                weights[i as usize][j as usize] = weight;
                weights[j as usize][i as usize] = weight;
            }
        }

        CompleteGraph {
            node_num: node_num,
            weights: weights,
        }
    }
}

fn sbm(
    graph: &CompleteGraph,
    alpha_0: f32,
    alpha_1: f32,
    time_step: f32,
    sub_step: u32,
    iteration: u32,
) -> Vec<f32> {
    let problem_size = graph.node_num as usize;
    let gamma_0 = 0.7 * alpha_0 / (graph.node_num as f32).sqrt();

    // コントロールパラメータ。0~1の範囲で変化する時間関数
    let mut control_alpha = 0.0;
    // すべてのオシレータに対するハミルトニアンにおける一般化座標
    let mut x: Vec<f32> = Vec::with_capacity(problem_size);
    for i in 0..graph.node_num {
        x[i as usize] = 0.0;
    }

    // すべてのオシレータに対するハミルトニアンにおける一般化運動量
    let mut p: Vec<f32> = Vec::with_capacity(problem_size);
    for i in 0..graph.node_num {
        let rng = rand::thread_rng();
        p[i as usize] = rng.gen_range(-0.1, 0.1);
    }

    let sub_time_step = time_step / (sub_step as f32);
    for l in 0..iteration {
        // オシレータの運動量の変化量を計算する。他のオシレータの相互作用を考慮する
        let mut delta_p: Vec<f32> = Vec::with_capacity(problem_size);
        for i in 0..problem_size {
            let mut sum = 0.0;
            for j in 0..problem_size {
                sum += graph.weights[i][j] * x[j];
            }
            let jx = time_step * gamma_0 * sum;
            delta_p.push(jx);
        }

        for i in 0..problem_size {
            p[i] += delta_p[i];
            for m in 0..sub_step {
                // 外部磁界の存在は無視する。外部磁界を考慮する場合、三項目が存在する。
                p[i] += sub_time_step * (-(alpha_0 - control_alpha) * x[i] - alpha_1 * (x[i] * *3));
                x[i] += sub_time_step * p[i]
            }
        }
    }
    // sgn(x_i)をとることで、イジングモデルにおけるスピンに変換する。つまり、+1 or -1
    // 通常の符号関数ではsgn(0) = 0であるが、初期値を0とした場合、解は更新されているはずなので、0である可能性を無視する。
    for i in 0..problem_size {
        if i < 0 {
            x[i] = -1.0;
        } else if i > 0 {
            x[i] = 1.0;
        } else {
            panic!(format!("Oscillator {} should be updated.", i));
        }
    }

    return x;
}

fn main() {
    println!("Hello, world")
}
