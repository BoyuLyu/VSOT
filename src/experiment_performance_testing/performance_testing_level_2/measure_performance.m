function [tamada_error_f1, our_method_error_f1, ofer_method_error_f1, dorkenwalk_error_f1,...
    tamada_error_precision, our_method_error_precision,ofer_method_error_precision,dorkenwalk_error_precision,...
        tamada_error_recall,our_method_error_recall,ofer_method_error_recall,dorkenwalk_error_recall,...
            tamada_error_iou,our_method_error_iou,ofer_method_error_iou, dorkenwalk_error_iou] = measure_performance(out1)
    tamada_error_f1 = out1.tamada_error_f1;
    tamada_error_precision = out1.tamada_error_precision;
    tamada_error_recall = out1.tamada_error_recall;
    tamada_error_iou = out1.tamada_error_iou;
    our_method_error_f1 = out1.our_method_error_f1;
    our_method_error_precision = out1.our_method_error_precision;
    our_method_error_recall = out1.our_method_error_recall;
    our_method_error_iou = out1.our_method_error_iou;
    ofer_method_error_f1 = out1.ofer_method_error_f1;
    ofer_method_error_precision = out1.ofer_method_error_precision;
    ofer_method_error_recall = out1.ofer_method_error_recall;
    ofer_method_error_iou = out1.ofer_method_error_iou;
    dorkenwalk_error_f1 = out1.dorkenwalk_error_f1;
    dorkenwalk_error_precision = out1.dorkenwalk_error_precision;
    dorkenwalk_error_recall = out1.dorkenwalk_error_recall;
    dorkenwalk_error_iou = out1.dorkenwalk_error_iou;
    
    nan_rm1 = find(isnan(tamada_error_f1)| isnan(our_method_error_f1) |isnan(ofer_method_error_f1)|isnan(dorkenwalk_error_f1));
    nan_rm2 = find(isnan(tamada_error_precision)| isnan(our_method_error_precision) |isnan(ofer_method_error_precision)|isnan(dorkenwalk_error_precision));
    nan_rm3 = find(isnan(tamada_error_recall)| isnan(our_method_error_recall) |isnan(ofer_method_error_recall)|isnan(dorkenwalk_error_recall));
    nan_rm4 = find(isnan(tamada_error_iou)| isnan(our_method_error_iou) |isnan(ofer_method_error_iou)|isnan(dorkenwalk_error_iou));
    
    tamada_error_f1(nan_rm1) = [];
    our_method_error_f1(nan_rm1) = [];
    ofer_method_error_f1(nan_rm1) = [];
    dorkenwalk_error_f1(nan_rm1) = [];
    tamada_error_precision(nan_rm2) = [];
    our_method_error_precision(nan_rm2) = [];
    ofer_method_error_precision(nan_rm2) = [];
    dorkenwalk_error_precision(nan_rm2) = [];
    tamada_error_recall(nan_rm3) = [];
    our_method_error_recall(nan_rm3) = [];
    ofer_method_error_recall(nan_rm3) = [];
    dorkenwalk_error_recall(nan_rm3) = [];
    tamada_error_iou(nan_rm4) = [];
    our_method_error_iou(nan_rm4) = [];
    ofer_method_error_iou(nan_rm4) = [];
    dorkenwalk_error_iou(nan_rm4) = [];
    
end