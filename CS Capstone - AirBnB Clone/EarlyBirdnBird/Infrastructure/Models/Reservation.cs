using Microsoft.AspNetCore.Mvc.ModelBinding.Validation;
using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    public class Reservation
    {
        //need to add foriegn key relationships
        [Key]
        public int Id { get; set; }



        [ForeignKey("RentalId")]
        public int Rental { get; set; }

        public double? IndividualTotal { get; set; }

        public double? SalesTax { get; set; }

        public double GrandTotal {  get; set; }

        public int? Guests { get; set; }

        public int Status { get; set; } // something like 1 = Pending, 2 = Approved, 3 = Cancelled ---  the intend is this can be used to filter 

        [DataType(DataType.Date)]
        public DateTime? BookedStart { get; set; }

        [DataType(DataType.Date)]
        public DateTime? BookedEnd { get; set; }

        [ForeignKey("AppUserId")]
        public AppUser? AppUser { get; set; }
    }
}
