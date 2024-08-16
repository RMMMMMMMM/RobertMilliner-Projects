using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Infrastructure.Models
{
    public class RentalBed
    {
        [Key]
        public int Id { get; set; }

        [ForeignKey("RentalId")]
        public int RentalId { get; set; }
        

        public Rental? Rental { get; set; }
        

        [ForeignKey("BedlId")] 
        public int BedId { get; set; }
        

        public Bed? Bed { get; set; }
    }
}