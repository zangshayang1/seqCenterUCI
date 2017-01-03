// JS

function jobSubmitFormCompleteCheck() {
    var clientName = document.forms["jobSubmitForm"]["clientName"].value;
    var clientEmail = document.forms["jobSubmitForm"]["clientEmail"].value;
    var selectService = document.forms["jobSubmitForm"]["selectService"].value;
    var sampleConcentration = document.forms["jobSubmitForm"]["sampleConcentration"].value;
    var selectSeqMode = document.forms["jobSubmitForm"]["selectSeqMode"].value;
    var selectReadsLength = document.forms["jobSubmitForm"]["selectReadsLength"].value;
    var targetDepth = document.forms["jobSubmitForm"]["targetDepth"].value;
    if (clientName.trim() == "" ||
        clientEmail.trim() == "" ||
        selectService.trim() == "" ||
        sampleConcentration.trim() == "" ||
        selectSeqMode.trim() == "" ||
        selectReadsLength.trim() == "" ||
        targetDepth.trim() == "") {
        
        alert("Failed. The form is incomplete.");
        return false;
    }
    return true;
}
